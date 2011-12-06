from basetest import *

from partfinder import alignment

class TestPhyml(PartitionFinderTestCase):

    def test_simple(self):
        test = """
5 10
spp1   actgactgaa
spp2   tctgtctgtt
spp3   agtgagtgaa
spp4   actaactaaa
spp5   acagacagaa
        """
        alignment.parse(test)

    def test_interleaved(self):
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
    unittest.main()
