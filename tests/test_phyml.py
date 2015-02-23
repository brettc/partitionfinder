from partfinder import alignment
from StringIO import StringIO

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

def test_simple():
    a = alignment.Alignment()
    a.parse(BASIC)
    assert a.species_count == 5
    assert a.sequence_length == 10
    assert a.data.shape == (5, 10)

    # write out and read back in
    out = StringIO()
    a.write_phylip(out)
    str_output = out.getvalue()
    del out

    b = alignment.Alignment()
    b.parse(str_output)

    assert (a.data == b.data).all()
    assert a.species == b.species

def test_interleaved():
    a = alignment.Alignment()
    a.parse(INTERLEAVED)
    assert a.species_count == 5
    assert a.sequence_length == 30
    assert a.data.shape == (5, 30)


class FakeSubset(object):
    def __init__(self, cols):
        self.columns = cols

def test_subset():
    a = alignment.Alignment()
    a.parse(INTERLEAVED)

    ss = FakeSubset([0, 1, 5, 7])
    b = alignment.SubsetAlignment(a, ss)
    assert b.sequence_length == len(ss.columns)
    for i, c in enumerate(ss.columns):
        assert (b.data[:, i] == a.data[:, c]).all()


