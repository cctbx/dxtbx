# LIBTBX_SET_DISPATCHER_NAME pytest
import sys

import pytest

# modify sys.argv so the command line help shows the right executable name
sys.argv[0] = "pytest"

sys.exit(pytest.main())
