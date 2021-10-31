import boost_adaptbx.boost.python

ext = boost_adaptbx.boost.python.import_ext("dxtbx_ext")

try:
    from .dxtbx_ext import *  # noqa: F401,F403,E402
except ModuleNotFoundError:
    from dxtbx_ext import *  # type: ignore # noqa: F401,F403,E402
