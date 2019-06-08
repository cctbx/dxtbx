from __future__ import absolute_import, division, print_function

import os
import sys
import xml.etree.ElementTree as ET

import dxtbx

tree = ET.parse("output.xml")
root = tree.getroot()
tcs = root.findall("testcase")

with open(os.path.join(dxtbx.__path__[0], ".python3-test-failures")) as fh:
    known_bad = set(testname.strip() for testname in fh.readlines() if testname.strip())
revised_bad = set()

acceptable_outcome = True

for test in tcs:
    test_name = test.attrib["file"] + "::" + test.attrib["name"]
    failure = test.find("failure") is not None or test.find("error") is not None
    if failure:
        if test_name in known_bad:
            revised_bad.add(test_name)
        else:
            print("Test %s has failed" % test_name)
            revised_bad.add(test_name)
            acceptable_outcome = False

for test in known_bad - revised_bad:
    print("Test %s did not fail. Please remove from known failure list" % test)

with open("revised-python3-test-failures", "w") as fh:
    fh.write("\n".join(sorted(revised_bad)))

sys.exit(0 if acceptable_outcome else 1)
