import os
from partfinder.raxml import likelihood_parser
import py.test

HERE = os.path.abspath(os.path.dirname(__file__))

raxml_likelihood_file = os.path.join(HERE, "RAxML_perSiteLLs.testing")

def test_likelihoods():
    likelihood_list = [[-2.705940], [-2.705940], [-2.698479], [-5.225005], \
        [-2.714398], [-2.698479], [-2.705940], [-2.705940], [-2.705940], \
        [-5.255713]]
    assert likelihood_parser(raxml_likelihood_file) == likelihood_list

def test_error():
    with py.test.raises(IOError):
        likelihood_parser("this_file_does_not_exist")
