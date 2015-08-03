import os
from partfinder import raxml
from partfinder.config import Configuration

HERE = os.path.abspath(os.path.dirname(__file__))
MISC_PATH = os.path.join(HERE, 'misc')


def test_parse_nucleotide():
    pth = os.path.join(MISC_PATH, 'raxml_nucleotide.output')
    c = Configuration().init(datatype='DNA', phylogeny_program="raxml")
    p = raxml.Parser(c)
    p.parse(open(pth).read())


def test_parse_aminoacid():
    pth = os.path.join(MISC_PATH, 'raxml_aminoacid.output')
    c = Configuration().init(datatype='protein', phylogeny_program="raxml")
    p = raxml.Parser(c)
    p.parse(open(pth).read())


