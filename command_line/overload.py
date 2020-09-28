import argparse

from dxtbx import load


def overload(image_file):
    i = load(image_file)
    data = i.get_raw_data()
    if not isinstance(data, tuple):
        data = (data,)
    detector = i.get_detector()
    for pid, (d, p) in enumerate(zip(data, detector)):
        if max(d) > p.get_trusted_range()[1]:
            return True
    return False


def run(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("image_file", nargs="+")
    options = parser.parse_args(args)
    for image_file in options.image_file:
        if overload(image_file):
            print(image_file)


if __name__ == "__main__":
    run()
