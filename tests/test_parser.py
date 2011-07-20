from basetest import *
from partfinder import config, parser
import os

class TestConfigFile(PartitionFinderTestCase):

    def test_config_creation(self):
        """Test loading a configuration file and checking values"""
        c = config.Configuration()
        c.load(os.path.join(self.cfg_path, 'test1.cfg'))

        self.assertEqual(c.alignment, 'test.phy')
        self.assertEqual(len(c.models), 56)

    def load_test(self, i):
        config_file = os.path.join(TestConfigFile.cfg_path, 'test%d.cfg' % i)
        c = config.Configuration()
        c.load(config_file)

# Dymanically add all separate files as tests
# Now we can just add 
# See here: http://stackoverflow.com/questions/1193909/pythons-unittest-and-dynamic-creation-of-test-cases
for i in range(1, 10):
    def ch(i):
        return lambda self: self.load_test(i)
    nm = "test_config_loading_%s" % i
    setattr(TestConfigFile, nm, ch(i))

if __name__ == '__main__':
    unittest.main()
