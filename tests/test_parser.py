from basetest import *
from partfinder import config, parser
import os


def load_config(f):
    config_file = os.path.join(CFG_PATH, f)
    c = config.Configuration()
    c.load(config_file)

def test_all_configs():
    cfg_files = os.listdir(CFG_PATH)
    cfg_files.sort()
    for f in cfg_files:
        yield load_config, f

if __name__ == '__main__':
    nose.runmodule()
