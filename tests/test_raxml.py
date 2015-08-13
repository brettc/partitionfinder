import os
from partfinder import raxml
from partfinder.config import Configuration
import numpy

HERE = os.path.abspath(os.path.dirname(__file__))
MISC_PATH = os.path.join(HERE, 'misc')

def test_parse_nucleotide():
    pth = os.path.join(MISC_PATH, 'raxml_nucleotide.output')
    c = Configuration().init(datatype='DNA', phylogeny_program="raxml")
    p = raxml.Parser(c)
    res = p.parse(open(pth).read())
    expected = numpy.array([
        0.315909, 
        0.232955,
        0.190909, 
        0.260227, 
    ])
    assert numpy.allclose(res.freqs[0], expected)

def test_parse_aminoacid():
    pth = os.path.join(MISC_PATH, 'raxml_aminoacid.output')
    c = Configuration().init(datatype='protein', phylogeny_program="raxml")
    p = raxml.Parser(c)
    res = p.parse(open(pth).read())
    expected = numpy.array([0.046704, 0.017460, 0.068093, 0.030991, 0.003492,
                            0.021388, 0.027499, 0.064164, 0.027499, 0.131820,
                            0.121781, 0.019642, 0.061545, 0.073330, 0.035356,
                            0.102139, 0.048014, 0.025753, 0.035356, 0.037975])
    assert numpy.allclose(res.freqs[0], expected)

def test_parse_lg4m():
    pth = os.path.join(MISC_PATH, 'raxml_aminoacid_LG4M+G.output')
    c = Configuration().init(datatype='protein', phylogeny_program="raxml")
    p = raxml.Parser(c)
    res = p.parse(open(pth).read())

    # Make sure we got the first block of rates
    expected = numpy.array([0.082276, 0.055172, 0.043853, 0.053484,
                            0.018957, 0.028152, 0.046679, 0.15781701,
                            0.033297, 0.028284, 0.054284, 0.025275,
                            0.023665, 0.041874, 0.063071, 0.066501,
                            0.065424, 0.023837, 0.038633, 0.049465])
    assert numpy.allclose(res.freqs[0], expected)
    
