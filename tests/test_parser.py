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
        # self.assertEqual(c.branchlengths, 'linked')
        # self.assertEqual(c.model_selection, 'bic')

    def test_config_file_loading(self):
        """Load all config files in 'cfg' folder"""
        for i in range(1, 10):
            config_file = os.path.join(self.cfg_path, 'test%d.cfg' % i)
            c = config.Configuration()
            c.load(config_file)

if __name__ == '__main__':
    unittest.main()
