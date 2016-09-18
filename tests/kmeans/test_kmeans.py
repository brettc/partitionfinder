import os
import pytest
from partfinder import main, util, analysis, config
from zipfile import ZipFile

HERE = os.path.abspath(os.path.dirname(__file__))

def test_grand():
    full_path = os.path.join(HERE, "Grande_2013")
    main.call_main("DNA", '--no-ml-tree --raxml "%s"' % full_path)
