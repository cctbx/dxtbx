import pickle

from dxtbx import IncorrectFormatError
from dxtbx.format.Format import Format


def test_pickle_incorrect_format():
    ex = IncorrectFormatError(Format, "Some message")
    pickle.loads(pickle.dumps(ex))
