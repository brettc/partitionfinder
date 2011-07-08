from basetest import *
from partfinder import config, parser
import os

class TestConfigFile(PartitionFinderTestCase):
    def test_config1(self):
        c = config.Configuration()
        p = parser.Parser(c)
        p.parse_file(os.path.join(self.cfg_path, 'test1.cfg'))

        self.assertEqual(c.alignment, 'test.phy')

    def test_config2(self):
        c = config.Configuration()
        p = parser.Parser(c)
        p.parse_file(os.path.join(self.cfg_path, 'test2.cfg'))

        # Could test for errors here too

if __name__ == '__main__':
    unittest.main()
