from basetest import *

from partfinder import alignment

def test_simple():
    test = """
5 10
spp1   actgactgaa
spp2   tctgtctgtt
spp3   agtgagtgaa
spp4   actaactaaa
spp5   acagacagaa
    """
    alignment.parse(test)

def test_interleaved():
    test = """
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
    alignment.parse(test)

if __name__ == '__main__':
    nose.runmodule()
