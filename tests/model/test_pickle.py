import six.moves.cPickle as pickle

from dxtbx.model import Beam, Detector, Goniometer, Panel, Scan


def pickle_then_unpickle(obj):
    """Pickle to a temp file then un-pickle."""
    return pickle.loads(pickle.dumps(obj))


def test_beam():
    """Test pickling the beam object."""
    obj1 = Beam((1, 1, 1))
    obj2 = pickle_then_unpickle(obj1)
    assert obj1 == obj2


def test_goniometer():
    """Test pickling the goniometer object."""
    obj1 = Goniometer()
    obj2 = pickle_then_unpickle(obj1)
    assert obj1 == obj2


def test_panel():
    """Test pickling the panel object."""
    obj1 = Panel()
    obj1.set_local_frame((1, 0, 0), (0, 1, 0), (0, 0, 1))
    obj2 = pickle_then_unpickle(obj1)
    assert obj1 == obj2


def test_detector():
    """Test pickling the detector object."""
    p = Panel()
    p.set_local_frame((1, 0, 0), (0, 1, 0), (0, 0, 1))
    obj1 = Detector(p)
    obj2 = pickle_then_unpickle(obj1)
    assert obj1 == obj2


def test_hierarchical_detector():
    """Test pickling the detector object."""
    p = Panel()
    p.set_local_frame((1, 0, 0), (0, 1, 0), (0, 0, 1))
    obj1 = Detector()
    root = obj1.hierarchy()
    root.add_panel(p)
    root.add_group()
    obj2 = pickle_then_unpickle(obj1)
    assert obj2.hierarchy()[0] == obj2[0]
    assert obj2.hierarchy()[0] in obj2
    assert obj2.hierarchy()[1].is_group()
    assert obj1 == obj2


def test_scan():
    """Test pickling the scan data object."""
    obj1 = Scan((1, 2), (1, 1))
    obj2 = pickle_then_unpickle(obj1)
    assert obj1 == obj2
