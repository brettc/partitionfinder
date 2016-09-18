import os
import pytest
from partfinder import main, util, analysis, config
from zipfile import ZipFile

HERE = os.path.abspath(os.path.dirname(__file__))

def test_DNA_entropy():
    full_path = os.path.join(HERE, "DNA_entropy")
    main.call_main("DNA", '--raxml --kmeans entropy --min-subset-size 10 "%s"' % full_path)

def test_DNA_tiger():
    full_path = os.path.join(HERE, "DNA_tiger")
    with pytest.raises(util.PartitionFinderError):
    	main.call_main("DNA", '--raxml --kmeans tiger --min-subset-size 10 "%s"' % full_path)

def test_morph_entropy():
    full_path = os.path.join(HERE, "morph_entropy")
    main.call_main("morphology", '--raxml --kmeans entropy --min-subset-size 1 "%s"' % full_path)

def test_morph_tiger():
    full_path = os.path.join(HERE, "morph_tiger")
    main.call_main("morphology", '--raxml --kmeans tiger --min-subset-size 1 "%s"' % full_path)
