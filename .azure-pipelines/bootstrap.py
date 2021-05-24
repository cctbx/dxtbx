import argparse
import multiprocessing.pool
import os
import re
import shutil
import socket as pysocket
import stat
import subprocess
import sys
import tarfile
import threading
import time
import zipfile
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

if sys.hexversion < 0x3060000:
    exit("Python 3.6+ required")

# Clean environment for subprocesses
clean_env = {
    key: value
    for key, value in os.environ.items()
    if key not in ("PYTHONPATH", "LD_LIBRARY_PATH", "DYLD_LIBRARY_PATH")
}

devnull = open(os.devnull, "wb")  # to redirect unwanted subprocess output
allowed_ssh_connections = {}
concurrent_git_connection_limit = threading.Semaphore(5)


def make_executable(filepath):
    if os.name == "posix":
        mode = os.stat(filepath).st_mode
        mode |= (mode & 0o444) >> 2  # copy R bits to X
        # r--r--r-- => 0o444
        os.chmod(filepath, mode)


def install_micromamba(python):
    """Download and install Micromamba"""
    if sys.platform.startswith("linux"):
        member = "bin/micromamba"
        url = "https://micromamba.snakepit.net/api/micromamba/linux-64/latest"
    elif sys.platform == "darwin":
        member = "bin/micromamba"
        url = "https://micromamba.snakepit.net/api/micromamba/osx-64/latest"
    elif os.name == "nt":
        member = "Library/bin/micromamba.exe"
        url = "https://micromamba.snakepit.net/api/micromamba/win-64/latest"
    else:
        raise NotImplementedError(f"Unsupported platform {os.name} / {sys.platform}")
    mamba_prefix = os.path.realpath("micromamba")
    clean_env["MAMBA_ROOT_PREFIX"] = mamba_prefix
    mamba = os.path.join(mamba_prefix, member.split("/")[-1])
    print(f"Downloading {url}:", end=" ")
    result = download_to_file(url, os.path.join(mamba_prefix, "micromamba.tar.bz2"))
    if result in (0, -1):
        sys.exit("Micromamba download failed")
    with tarfile.open(
        os.path.join(mamba_prefix, "micromamba.tar.bz2"), "r:bz2"
    ) as tar, open(mamba, "wb") as fh:
        fh.write(tar.extractfile(member).read())
    make_executable(mamba)

    # verify micromamba works and check version
    conda_info = subprocess.check_output([mamba, "--version"], env=clean_env)
    if sys.version_info.major > 2:
        conda_info = conda_info.decode("latin-1")
    print("Using Micromamba version", conda_info.strip())

    # identify packages required for environment
    filename = os.path.join("modules", "dxtbx", ".azure-pipelines", "ci-conda-env.txt")
    if not os.path.isfile(filename):
        raise RuntimeError(f"The environment file {filename} is not available")

    # install a new environment or update an existing one
    prefix = os.path.realpath("conda_base")
    if os.path.exists(prefix):
        command = "install"
        text_messages = ["Updating", "update of"]
    else:
        command = "create"
        text_messages = ["Installing", "installation into"]

    command_list = [
        mamba,
        "--no-env",
        "--no-rc",
        "--prefix",
        prefix,
        "--root-prefix",
        mamba_prefix,
        command,
        "--file",
        filename,
        "--yes",
        "--channel",
        "conda-forge",
        "--override-channels",
        "python=%s" % python,
    ]

    print(
        "{text} dials environment from {filename} with Python {python}".format(
            text=text_messages[0], filename=filename, python=python
        )
    )
    for retry in range(5):
        retry += 1
        try:
            run_command(
                command=command_list,
                workdir=".",
            )
        except Exception:
            print(
                """
*******************************************************************************
There was a failure in constructing the conda environment.
Attempt {retry} of 5 will start {retry} minute(s) from {t}.
*******************************************************************************
""".format(
                    retry=retry, t=time.asctime()
                )
            )
            time.sleep(retry * 60)
        else:
            break
    else:
        sys.exit(
            """
The conda environment could not be constructed. Please check that there is a
working network connection for downloading conda packages.
"""
        )
    print("Completed {text}:\n  {prefix}".format(text=text_messages[1], prefix=prefix))
    with open(os.path.join(prefix, ".condarc"), "w") as fh:
        fh.write(
            """
changeps1: False
channels:
  - conda-forge
""".lstrip()
        )


def run_command(command, workdir):
    print("Running {} (in {})".format(" ".join(command), workdir))
    workdir = os.path.abspath(workdir)
    try:
        os.makedirs(workdir)
    except OSError:
        pass
    try:
        p = subprocess.Popen(args=command, cwd=workdir, env=clean_env)
    except Exception as e:
        if isinstance(e, OSError):
            if e.errno == 2:
                executable = os.path.normpath(os.path.join(workdir, command[0]))
                raise RuntimeError("Could not run %s: File not found" % executable)
        if "child_traceback" in dir(e):
            print("Calling subprocess resulted in error; ", e.child_traceback)
        raise e

    try:
        p.wait()
    except KeyboardInterrupt:
        print("\nReceived CTRL+C, trying to stop subprocess...\n")
        p.terminate()
        raise
    if p.returncode:
        sys.exit("Process failed with return code %s" % p.returncode)


def run_indirect_command(command, args):
    print("(via conda environment) " + command)
    if os.name == "nt":
        filename = os.path.join("build", "indirection.cmd")
        with open(filename, "w") as fh:
            fh.write("call %s\\conda_base\\condabin\\activate.bat\r\n" % os.getcwd())
            fh.write("shift\r\n")
            fh.write("%*\r\n")
        command = command + ".bat"
        indirection = ["cmd.exe", "/C", "indirection.cmd"]
    else:
        filename = os.path.join("build", "indirection.sh")
        with open(filename, "w") as fh:
            fh.write("#!/bin/bash\n")
            fh.write("source %s/conda_base/etc/profile.d/conda.sh\n" % os.getcwd())
            fh.write("conda activate %s/conda_base\n" % os.getcwd())
            fh.write('"$@"\n')
        make_executable(filename)
        indirection = ["./indirection.sh"]
    run_command(
        command=indirection + [command] + args,
        workdir="build",
    )


def download_to_file(url, file, quiet=False):
    """Downloads a URL to file. Returns the file size.
    Returns -1 if the downloaded file size does not match the expected file
    size
    Returns -2 if the download is skipped due to the file at the URL not
    being newer than the local copy (identified by matching timestamp and
    size)
    """

    # Create directory structure if necessary
    if os.path.dirname(file):
        try:
            os.makedirs(os.path.dirname(file))
        except Exception:
            pass

    localcopy = os.path.isfile(file)

    try:
        from ssl import SSLError
    except ImportError:
        SSLError = None

    # Open connection to remote server
    try:
        url_request = Request(url)
        if localcopy:
            # Shorten timeout to 7 seconds if a copy of the file is already present
            socket = urlopen(url_request, None, 7)
        else:
            socket = urlopen(url_request)
    except SSLError as e:
        # This could be a timeout
        if localcopy:
            # Download failed for some reason, but a valid local copy of
            # the file exists, so use that one instead.
            if not quiet:
                print(str(e))
            return -2
        # otherwise pass on the error message
        raise
    except (pysocket.timeout, HTTPError) as e:
        if localcopy:
            # Download failed for some reason, but a valid local copy of
            # the file exists, so use that one instead.
            if not quiet:
                print(str(e))
            return -2
        # otherwise pass on the error message
        raise
    except URLError as e:
        if localcopy:
            # Download failed for some reason, but a valid local copy of
            # the file exists, so use that one instead.
            if not quiet:
                print(str(e))
            return -2
        # if url fails to open, try using curl
        # temporary fix for old OpenSSL in system Python on macOS
        # https://github.com/cctbx/cctbx_project/issues/33
        command = ["/usr/bin/curl", "--http1.0", "-fLo", file, "--retry", "5", url]
        subprocess.call(command)
        socket = None  # prevent later socket code from being run
        try:
            received = os.path.getsize(file)
        except OSError:
            raise RuntimeError("Download failed")

    if socket is not None:
        try:
            file_size = int(socket.info().get("Content-Length"))
        except Exception:
            file_size = 0

        remote_mtime = 0
        try:
            remote_mtime = time.mktime(socket.info().getdate("last-modified"))
        except Exception:
            pass

        if file_size > 0:
            if remote_mtime > 0:
                # check if existing file matches remote size and timestamp
                try:
                    (
                        mode,
                        ino,
                        dev,
                        nlink,
                        uid,
                        gid,
                        size,
                        atime,
                        mtime,
                        ctime,
                    ) = os.stat(file)
                    if (size == file_size) and (remote_mtime == mtime):
                        if not quiet:
                            print("local copy is current")
                        socket.close()
                        return -2
                except Exception:
                    # proceed with download if timestamp/size check fails for any reason
                    pass

            hr_size = (file_size, "B")
            if hr_size[0] > 500:
                hr_size = (hr_size[0] / 1024, "kB")
            if hr_size[0] > 500:
                hr_size = (hr_size[0] / 1024, "MB")
            if not quiet:
                print("%.1f %s" % hr_size)
                print("    [0%", end="")
                sys.stdout.flush()  # becomes print(flush=True) when we move to 3.3+

        received = 0
        block_size = 8192
        progress = 1
        # Write to the file immediately so we can empty the buffer
        tmpfile = file + ".tmp"

        with open(tmpfile, "wb") as fh:
            while True:
                block = socket.read(block_size)
                received += len(block)
                fh.write(block)
                if file_size > 0 and not quiet:
                    while (100 * received / file_size) > progress:
                        progress += 1
                        if (progress % 20) == 0:
                            print(progress, end="%")
                            sys.stdout.flush()  # becomes print(flush=True) when we move to 3.3+
                        elif (progress % 2) == 0:
                            print(".", end="")
                            sys.stdout.flush()  # becomes print(flush=True) when we move to 3.3+
                if not block:
                    break
        socket.close()

        if not quiet:
            if file_size > 0:
                print("]")
            else:
                print("%d kB" % (received / 1024))
            sys.stdout.flush()  # becomes print(flush=True) when we move to 3.3+

        # Do not overwrite file during the download. If a download temporarily fails we
        # may still have a clean, working (yet older) copy of the file.
        shutil.move(tmpfile, file)

        if (file_size > 0) and (file_size != received):
            return -1

        if remote_mtime > 0:
            # set file timestamp if timestamp information is available
            st = os.stat(file)
            atime = st[stat.ST_ATIME]  # current access time
            os.utime(file, (atime, remote_mtime))

    return received


def unzip(archive, directory, trim_directory=0):
    """unzip a file into a directory."""
    if not zipfile.is_zipfile(archive):
        raise Exception(
            f"Cannot install {directory}: {archive} is not a valid .zip file"
        )
    z = zipfile.ZipFile(archive, "r")
    for member in z.infolist():
        is_directory = member.filename.endswith("/")
        filename = os.path.join(*member.filename.split("/")[trim_directory:])
        if filename != "":
            filename = os.path.normpath(filename)
            if "../" in filename:
                raise Exception(
                    f"Archive {archive} contains invalid filename {filename}"
                )
            filename = os.path.join(directory, filename)
            upperdirs = os.path.dirname(filename)
            try:
                if is_directory and not os.path.exists(filename):
                    os.makedirs(filename)
                elif upperdirs and not os.path.exists(upperdirs):
                    os.makedirs(upperdirs)
            except Exception:
                pass
            if not is_directory:
                with z.open(member) as source:
                    with open(filename, "wb") as target:
                        shutil.copyfileobj(source, target)

                # Preserve executable permission, if set
                unix_executable = member.external_attr >> 16 & 0o111
                # rwxrwxrwx => --x--x--x => 0o111
                if unix_executable:
                    make_executable(filename)
    z.close()


def set_git_repository_config_to_rebase(config):
    with open(config) as fh:
        cfg = fh.readlines()

    branch, remote, rebase = False, False, False
    insertions = []
    for n, line in enumerate(cfg):
        if line.startswith("["):
            if branch and remote and not rebase:
                insertions.insert(0, (n, branch))
            if line.startswith("[branch"):
                branch = line.split('"')[1]
            else:
                branch = False
            remote, rebase = False, False
        if re.match(r"remote\s*=", line.strip()):
            remote = True
        if re.match(r"rebase\s*=", line.strip()):
            rebase = True
    if branch and remote and not rebase:
        insertions.insert(0, (n + 1, branch))
    for n, branch in insertions:
        cfg.insert(n, "\trebase = true\n")
    with open(config, "w") as fh:
        fh.write("".join(cfg))


def git(module, git_available, ssh_available, settings):
    """Retrieve a git repository, either by running git directly
    or by downloading and unpacking an archive.
    """
    destination = os.path.join("modules", module)

    if os.path.exists(destination):
        if not os.path.exists(os.path.join(destination, ".git")):
            return module, "WARNING", "Existing non-git directory -- skipping"
        if not git_available:
            return module, "WARNING", "Cannot update module, git command not found"

        with open(os.path.join(destination, ".git", "HEAD")) as fh:
            if fh.read(4) != "ref:":
                return (
                    module,
                    "WARNING",
                    "Cannot update existing git repository! You are not on a branch.\n"
                    "This may be legitimate when run eg. via Jenkins, but be aware that you cannot commit any changes",
                )

        with concurrent_git_connection_limit:
            p = subprocess.Popen(
                args=["git", "pull", "--rebase"],
                cwd=destination,
                env=clean_env,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )
            # This may fail for unclean trees and merge problems. In this case manual
            # user intervention will be required.
            # For the record, you can clean up the tree and *discard ALL changes* with
            #   git reset --hard origin/master
            #   git clean -dffx
            try:
                output, _ = p.communicate()
                output = output.decode("latin-1")
            except KeyboardInterrupt:
                print("\nReceived CTRL+C, trying to terminate subprocess...\n")
                p.terminate()
                raise
        if p.returncode:
            return (
                module,
                "WARNING",
                "Cannot update existing git repository! Unclean tree or merge problems.\n"
                + output,
            )
        # Show the hash for the checked out commit for debugging purposes
        p = subprocess.Popen(
            args=["git", "rev-parse", "HEAD", "--abbrev-ref", "HEAD"],
            cwd=destination,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        output, _ = p.communicate()
        output = output.decode("latin-1")
        if p.returncode:
            return module, "WARNING", "Cannot get git repository revision\n" + output
        output = output.split()
        if len(output) == 2:
            return (
                module,
                "OK",
                "Checked out revision {} ({})".format(output[0], output[1]),
            )
        return module, "OK", "Checked out revision " + output[0].strip()

    try:
        os.makedirs("modules")
    except OSError:
        pass

    remote_branch = settings.get("branch-remote", settings["branch-local"])

    if not git_available:
        # Fall back to downloading a static archive
        url = "https://github.com/{}/archive/{}.zip".format(
            settings.get("effective-repository", settings.get("base-repository")),
            remote_branch,
        )
        filename = os.path.join("modules", f"{module}-{remote_branch}.zip")
        try:
            download_to_file(url, filename, quiet=True)
        except Exception:
            print("Error downloading", url)
            raise
        unzip(filename, destination, trim_directory=1)
        return module, "OK", "Downloaded branch %s from static archive" % remote_branch

    if ssh_available:
        remote_pattern = "git@github.com:%s.git"
    else:
        remote_pattern = "https://github.com/%s.git"

    secondary_remote = settings.get("effective-repository") and (
        settings["effective-repository"] != settings.get("base-repository")
    )
    direct_branch_checkout = []
    if not secondary_remote and remote_branch == settings["branch-local"]:
        direct_branch_checkout = ["-b", remote_branch]

    with concurrent_git_connection_limit:
        p = subprocess.Popen(
            args=["git", "clone", "--recursive"]
            + direct_branch_checkout
            + [
                remote_pattern
                % settings.get("base-repository", settings.get("effective-repository")),
                module,
            ],
            cwd="modules",
            env=clean_env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        try:
            output, _ = p.communicate()
            output = output.decode("latin-1")
        except KeyboardInterrupt:
            print("\nReceived CTRL+C, trying to terminate subprocess...\n")
            p.terminate()
            raise
        if p.returncode:
            return (module, "ERROR", "Cannot checkout git repository\n" + output)

    if secondary_remote:
        returncode = subprocess.call(
            [
                "git",
                "remote",
                "add",
                "upstream",
                remote_pattern % settings["effective-repository"],
            ],
            cwd=destination,
            env=clean_env,
            stdout=devnull,
            stderr=devnull,
        )
        if returncode:
            return (
                module,
                "ERROR",
                "Could not add upstream remote to repository. Repository may be in invalid state!",
            )
        with concurrent_git_connection_limit:
            returncode = subprocess.call(
                ["git", "fetch", "upstream"],
                cwd=destination,
                env=clean_env,
                stdout=devnull,
                stderr=devnull,
            )
        if returncode:
            return (
                module,
                "ERROR",
                "Could not fetch upstream repository %s. Repository may be in invalid state!"
                % settings["effective-repository"],
            )

    set_git_repository_config_to_rebase(os.path.join(destination, ".git", "config"))

    if not direct_branch_checkout:
        # set up the local branch with tracking
        returncode = subprocess.call(
            [
                "git",
                "checkout",
                "-B",
                settings["branch-local"],
                "--track",
                "{}/{}".format(
                    "upstream" if secondary_remote else "origin", remote_branch
                ),
            ],
            cwd=destination,
            env=clean_env,
            stdout=devnull,
            stderr=devnull,
        )
        if returncode:
            return (
                module,
                "ERROR",
                "Could not check out alternate branch %s. Repository may be in invalid state!"
                % remote_branch,
            )

    # Show the hash for the checked out commit for debugging purposes
    p = subprocess.Popen(
        args=["git", "rev-parse", "HEAD"],
        cwd=destination,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    output, _ = p.communicate()
    output = output.decode("latin-1")
    if p.returncode:
        return (
            module,
            "WARNING",
            "Cannot get git repository revision\n" + output,
        )
    git_status = settings["branch-local"]
    if settings["branch-local"] != remote_branch:
        git_status += " tracking " + remote_branch
    if secondary_remote:
        git_status += " at " + settings["effective-repository"]
    return module, "OK", f"Checked out revision {output.strip()} ({git_status})"


def update_sources(options):
    try:
        git_available = not subprocess.call(
            ["git", "--version"], stdout=devnull, stderr=devnull
        )
    except OSError:
        git_available = False
    ssh_available = False
    if git_available:
        try:
            returncode = subprocess.call(
                [
                    "ssh",
                    "-oBatchMode=yes",
                    "-oStrictHostKeyChecking=no",
                    "-T",
                    "git@github.com",
                ],
                stdout=devnull,
                stderr=devnull,
            )
            # SSH errors lead to 255
            ssh_available = returncode in (0, 1)
        except OSError:
            pass

    repositories = {
        source.split("/")[1]: {"base-repository": source, "branch-local": "master"}
        for source in (
            "cctbx/cctbx_project",
            "dials/cbflib",
        )
    }
    repositories["cctbx_project"] = {
        "base-repository": "cctbx/cctbx_project",
        "effective-repository": "dials/cctbx",
        "branch-remote": "master",
        "branch-local": "stable",
    }
    repositories["dxtbx"] = {
        "base-repository": "cctbx/dxtbx",
        "branch-local": "main",
    }

    for source, setting in options.branch:
        if source not in repositories:
            sys.exit("Unknown repository %s" % source)
        setting = re.match(
            r"^(?:(\w+/\w+)@?)?([a-zA-Z0-9._\-]+)?(?::([a-zA-Z0-9._\-]+))?$", setting
        )
        if not setting:
            sys.exit("Could not parse the branch setting for repository %s" % source)
        _repository, _branch_remote, _branch_local = setting.groups()
        if _repository:
            repositories[source] = {
                "base-repository": _repository,
                "branch-remote": _branch_remote or "master",
                "branch-local": _branch_local or _branch_remote or "master",
            }
        elif _branch_remote:
            repositories[source]["branch-remote"] = _branch_remote
            repositories[source]["branch-local"] = _branch_local or _branch_remote
        elif _branch_local:
            repositories[source]["branch-local"] = _branch_local

    def _git_fn(repository):
        return git(repository, git_available, ssh_available, repositories[repository])

    success = True
    update_pool = multiprocessing.pool.ThreadPool(20)
    try:
        for result in update_pool.imap_unordered(_git_fn, repositories):
            module, result, output = result
            output = (result + " - " + output).replace(
                "\n", "\n" + " " * (len(module + result) + 5)
            )
            if os.name == "posix" and sys.stdout.isatty():
                if result == "OK":
                    output = "\x1b[32m" + output + "\x1b[0m"
                elif result == "WARNING":
                    output = "\x1b[33m" + output + "\x1b[0m"
                elif result == "ERROR":
                    output = "\x1b[31m" + output + "\x1b[0m"
                    success = False
            print(module + ": " + output)
    except KeyboardInterrupt:
        update_pool.terminate()
        sys.exit("\naborted with Ctrl+C")
    except Exception:
        update_pool.terminate()
        raise
    update_pool.close()
    update_pool.join()
    if not success:
        sys.exit("\nFailed to update one or more repositories")


def refresh_build():
    print("Running libtbx.refresh")
    dispatch_extension = ".bat" if os.name == "nt" else ""
    run_command(
        [os.path.join("build", "bin", "libtbx.refresh" + dispatch_extension)],
        workdir=".",
    )


def configure_build(options):
    if os.name == "nt":
        conda_python = os.path.join(os.getcwd(), "conda_base", "python.exe")
    else:
        conda_python = os.path.join("..", "conda_base", "bin", "python")

    with open("dials", "w") as fh:
        # ensure we write a new-style environment setup script
        fh.write(
            """
#!/bin/bash

if [ -n "$BASH_SOURCE" ]; then
  LIBTBX_BUILD="$(dirname -- "${BASH_SOURCE[0]}")/build"
else
  LIBTBX_BUILD="%s"
fi

# make path absolute and resolve symlinks
LIBTBX_BUILD=$(cd -P -- "${LIBTBX_BUILD}" && pwd -P)

# enable conda environment
source ${LIBTBX_BUILD}/../conda_base/etc/profile.d/conda.sh
conda activate $(dirname -- "${LIBTBX_BUILD}")/conda_base

# prepend cctbx /build/bin directory to PATH
PATH="${LIBTBX_BUILD}/bin:${PATH}"
export PATH

unset LIBTBX_BUILD
""".lstrip()
            % os.path.join(os.getcwd(), "build")
        )
    print("Setting up build directory")
    run_command(
        command=[
            conda_python,
            os.path.join("..", "modules", "cctbx_project", "libtbx", "configure.py"),
            "--skip_phenix_dispatchers",
            "cctbx",
            "cbflib",
            "dxtbx",
            "scitbx",
            "libtbx",
            "iotbx",
            "--compiler=conda",
            "--enable_cxx11",
            "--use_conda",
        ],
        workdir="build",
    )


def make_build():
    try:
        nproc = len(os.sched_getaffinity(0))
    except AttributeError:
        nproc = multiprocessing.cpu_count()
    run_indirect_command(os.path.join("bin", "libtbx.scons"), args=["-j", str(nproc)])
    # run build again to make sure everything is built
    run_indirect_command(os.path.join("bin", "libtbx.scons"), args=["-j", str(nproc)])


def repository_at_tag(string):
    try:
        repository, tag = string.split("@", 1)
        return (repository, tag)
    except ValueError:
        raise argparse.ArgumentTypeError(
            "%s does not follow the repository@branch format" % string
        )


class Choices(tuple):
    # Python bug https://bugs.python.org/issue27227, https://bugs.python.org/issue9625
    def __new__(cls, *args, **kwargs):
        x = tuple.__new__(cls, *args, **kwargs)
        Choices.__init__(x, *args, **kwargs)
        return x

    def __init__(self, *args, **kwargs):
        self.default = []

    def __contains__(self, item):
        return tuple.__contains__(self, item) or item is self.default


def run():
    description = "A bootstrap helper for dxtbx CI builds"

    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    action_choices = Choices(("update", "base", "build"))
    action_choices.default = ["update", "base", "build"]
    parser.add_argument(
        "actions",
        nargs="*",
        help="actions for building",
        choices=action_choices,
        default=action_choices.default,
    )
    parser.add_argument(
        "--python",
        help="Install this minor version of Python (default: %(default)s)",
        default="3.8",
        choices=("3.6", "3.7", "3.8", "3.9"),
    )
    parser.add_argument(
        "--branch",
        type=repository_at_tag,
        action="append",
        default=[],
        help=(
            "during 'update' step when a repository is newly cloned set it to a given branch."
            "Specify as repository@branch, eg. 'dials@dials-next'"
        ),
    )
    parser.add_argument(
        "--clean",
        help="Remove temporary conda environments and package caches after installation",
        default=False,
        action="store_true",
    )
    options = parser.parse_args()

    # Add sources
    if "update" in options.actions:
        update_sources(options)

    # Build base packages
    if "base" in options.actions:
        install_micromamba(options.python)
        if options.clean:
            shutil.rmtree(os.path.realpath("micromamba"))

    # Configure, make
    if "build" in options.actions:
        configure_build(options)
        make_build()
        refresh_build()

    print("\nBootstrap success: %s" % ", ".join(options.actions))


if __name__ == "__main__":
    run()
