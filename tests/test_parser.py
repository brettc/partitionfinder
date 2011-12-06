from basetest import *
from partfinder import config, parser
import os

class TestConfigFile(PartitionFinderTestCase):

    def test_config_creation(self):
        """Test loading a configuration file and checking values"""
        c = config.Configuration()
        c.load(os.path.join(CFG_PATH, 'test1.cfg'))
        print c.alignment
        print c.user_tree_topology

        self.assertEqual(c.alignment, 'test.phy')
        self.assertEqual(len(c.models), 56)

    def load_test(self, f):
        config_file = os.path.join(CFG_PATH, f)
        c = config.Configuration()
        c.load(config_file)

# Dymanically add all separate files as tests
# Now we can just add new files
# See here: http://stackoverflow.com/questions/1193909/pythons-unittest-and-dynamic-creation-of-test-cases
# cfg_files = os.listdir(CFG_PATH)
# cfg_files.sort()
# for f in cfg_files:
    # def ch(f):
        # return lambda self: self.load_test(f)
    # setattr(TestConfigFile, f, ch(f))

if __name__ == '__main__':
    unittest.main()
