from basetest import *
from partfinder import config, parser
import os

class TestConfigFile(PartitionFinderTestCase):

    def test_config_creation(self):
        c = config.Configuration()
        p = parser.Parser(c)
        p.parse_file(os.path.join(self.cfg_path, 'test1.cfg'))

        self.assertEqual(c.alignment, 'test.phy')
        self.assertEqual(len(c.models), 56)
        self.assertEqual(c.branchlengths, 'linked')
        self.assertEqual(c.model_selection, 'bic')

    def test_config_file_loading(self):
        '''All of these config files should load without errors'''
        for i in range(1,10):
            c = config.Configuration()
            p = parser.Parser(c)
            config_file = 'test%d.cfg' %(i)
            print "Loading:", config_file
            p.parse_file(os.path.join(self.cfg_path, config_file))


        # Could test for errors here too

if __name__ == '__main__':
    unittest.main()
