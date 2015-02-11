# from partfinder.phyml import likelihood_parser
# from math import log
import os
# import py.test
import numpy as np

HERE = os.path.abspath(os.path.dirname(__file__))

phyml_file_4 = os.path.join(HERE, "free_rates_4_cat.phy_phyml_lk.txt")
phyml_file_8 = os.path.join(HERE, "gamma_8_cat.phy_phyml_lk.txt")
phyml_file_likelihoods = os.path.join(HERE, "simulated_TRUE.phy_phyml_lk.txt")

def test_something():
    x = np.genfromtxt(phyml_file_8, skiprows=7, usecols=range(1, 8), dtype=float)
    print x.shape
    print x[0]
    print x.T

if __name__=='__main__':
    test_something()
