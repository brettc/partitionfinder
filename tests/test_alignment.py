import pytest
import fnmatch
import os
from StringIO import StringIO
from partfinder.alignment import Alignment, SubsetAlignment, AlignmentError

HERE = os.path.abspath(os.path.dirname(__file__))
MISC_PATH = os.path.join(HERE, 'misc')

BASIC = """
5 10
spp1   actgactgaa
spp2   tctgtctgtt
spp3   agtgagtgaa
spp4   actaactaaa
spp5   acagacagaa
"""

INTERLEAVED = """
5 30
spp1   actgactgaa
spp2   tctgtctgtt
spp3   agtgagtgaa
spp4   actaactaaa
spp5   acagacagaa

actgactgaa
tctgtctgtt
agtgagtgaa
actaactaaa
acagacagaa

actgactgaa
tctgtctgtt
agtgagtgaa
actaactaaa
acagacagaa
"""

TOO_FEW_SPECIES = """
3 2
a AT
b AT
c AT
d AC
"""

TOO_MANY_SPECIES = """
5 2
a AT
b AT
c AT
d AC
"""


def write_and_get_stream(align):
    # write out and read back in
    out = StringIO()
    align.write_phylip(out)
    str_output = out.getvalue()
    del out
    return str_output


def test_simple():
    a = Alignment()
    a.parse(BASIC)
    assert a.species_count == 5
    assert a.sequence_length == 10
    assert a.data.shape == (5, 10)

    b = Alignment()
    b.parse(write_and_get_stream(a))

    assert (a.data == b.data).all()
    assert a.species == b.species


def test_interleaved():
    a = Alignment()
    a.parse(INTERLEAVED)
    assert a.species_count == 5
    assert a.sequence_length == 30
    assert a.data.shape == (5, 30)


def test_too_few_species(caplog):
    a = Alignment()
    with pytest.raises(AlignmentError):
        a.parse(TOO_FEW_SPECIES)
    assert "too many species" in caplog.text()


def test_too_many_species(caplog):
    a = Alignment()
    with pytest.raises(AlignmentError):
        a.parse(TOO_MANY_SPECIES)
    assert "Phyml format error" in caplog.text()


class FakeSubset(object):
    def __init__(self, cols):
        self.columns = cols


def test_subset():
    a = Alignment()
    a.parse(INTERLEAVED)

    ss = FakeSubset([0, 1, 5, 7])
    b = SubsetAlignment(a, ss)
    assert b.sequence_length == len(ss.columns)
    for i, c in enumerate(ss.columns):
        assert (b.data[:, i] == a.data[:, c]).all()


def generate_phyml_paths():
    paths = []
    for pth in os.listdir(MISC_PATH):
        if fnmatch.fnmatch(pth, "*.phy"):
            paths.append(os.path.join(MISC_PATH, pth))

    return paths


@pytest.fixture(params=generate_phyml_paths())
def phyml_path(request):
    return request.param


def test_load_and_save_phy(phyml_path):
    a = Alignment()
    a.read(phyml_path)
    b = Alignment()
    b.parse(write_and_get_stream(a))

    assert (a.data == b.data).all()
    assert a.species == b.species
