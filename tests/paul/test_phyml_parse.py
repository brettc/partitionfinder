from partfinder.phyml import likelihood_parser
from math import log
import os
import py.test

HERE = os.path.abspath(os.path.dirname(__file__))

phyml_file_4 = os.path.join(HERE, "free_rates_4_cat.phy_phyml_lk.txt")
phyml_file_8 = os.path.join(HERE, "gamma_8_cat.phy_phyml_lk.txt")
phyml_file_likelihoods = os.path.join(HERE, "simulated_TRUE.phy_phyml_lk.txt")

def test_4_cat():
    likelihood_tuple_4 = ([[log(2.07027e-12)], [log(1.8652e-07)], \
        [log(4.48873e-15)], [log(3.38958e-10)], [log(8.29969e-17)], [log(9.24579e-09)], [log(3.43996e-10)], \
        [log(4.43262e-13)], [log(3.42513e-11)], [log(1.15506e-11)]], \
        [[log(1.3895e-19), log(6.2676e-12), log(1.2534e-12), log(1.21786e-15)], \
        [log(2.05811e-19), log(6.73481e-07), log(4.14575e-09), log(7.97623e-14)], \
        [log(1.37274e-19), log(7.11221e-15), log(9.11826e-15), log(9.21848e-17)], \
        [log(1.31413e-19), log(1.18939e-09), log(4.20659e-11), log(5.86537e-15)], \
        [log(1.11587e-19), log(3.1672e-17), log(2.52183e-16), log(1.9722e-17)], \
        [log(1.59891e-19), log(3.31101e-08), log(4.79946e-10), log(2.59524e-14)], \
        [log(2.1917e-19), log(1.19544e-09), log(5.43128e-11), log(1.22969e-14)], \
        [log(1.1447e-19), log(1.32148e-12), log(2.8874e-13), log(3.7386e-16)], \
        [log(1.70149e-19), log(1.14227e-10), log(1.02103e-11), log(4.05239e-15)], \
        [log(1.28107e-19), log(3.86378e-11), log(3.32642e-12), log(1.46151e-15)]])
    assert likelihood_parser(phyml_file_4) == likelihood_tuple_4

def test_8_cat():
    likelihood_tuple_8 = ([[log(1.36373e-12)], [log(2.87895e-07)], [log(3.67318e-15)], \
        [log(2.75073e-10)], [log(8.1444e-17)], [log(9.90893e-09)], [log(2.54015e-10)], \
        [log(3.07761e-13)], [log(2.23513e-11)], [log(7.8803e-12)]], \
        [[log(3.60785e-12), log(5.44944e-12), log(1.58289e-12), log(2.45474e-13), \
        log(2.29162e-14), log(1.22494e-15), log(3.35077e-17), log(4.90601e-19)], \
        [log(2.16944e-06), log(1.27515e-07), log(5.98776e-09), log(2.18854e-10), \
        log(5.42051e-12), log(7.49452e-14), log(4.51885e-16), log(1.18074e-18)], \
        [log(1.63466e-15), log(1.32314e-14), log(1.01795e-14), log(3.5268e-15), \
        log(7.16337e-16), log(8.94944e-17), log(6.8864e-18), log(3.49039e-19)], \
        [log(1.64216e-09), log(4.95499e-10), log(5.82998e-11), log(4.41172e-12), \
        log(2.1336e-13), log(5.88679e-15), log(7.95885e-17), log(5.64329e-19)], \
        [log(2.82279e-18), log(1.22459e-16), log(2.45719e-16), log(1.83404e-16), \
        log(7.52518e-17), log(1.87904e-17), log(2.85338e-18), log(2.52601e-19)], \
        [log(6.94317e-08), log(9.12817e-09), log(6.74745e-10), log(3.55841e-11), \
        log(1.25053e-12), log(2.52239e-14), log(2.34373e-16), log(8.81134e-19)], \
        [log(1.44787e-09), log(5.08831e-10), log(6.89929e-11), log(6.06564e-12), \
        log(3.4634e-13), log(1.14327e-14), log(1.77155e-16), log(1.07488e-18)], \
        [log(7.83485e-13), log(1.23463e-12), log(3.7553e-13), log(6.17754e-14), \
        log(6.26968e-15), log(3.81829e-16), log(1.31831e-17), log(3.20361e-19)], \
        [log(9.637e-11), log(6.79003e-11), log(1.29863e-11), log(1.44982e-12), \
        log(1.00113e-13), log(3.93728e-15), log(7.65775e-17), log(7.18937e-19)], \
        [log(3.47616e-11), log(2.33762e-11), log(4.37891e-12), log(4.89457e-13), \
        log(3.47469e-14), log(1.45604e-15), log(3.20495e-17), log(4.30831e-19)]])
    assert likelihood_parser(phyml_file_8) == likelihood_tuple_8

def test_1_cat():
    likelihood_list = [[log(0.0565097)], [log(0.0565097)], [log(0.0580009)], \
        [log(0.00728218)], [log(0.0569829)], [log(0.0580009)], [log(0.0565097)], \
        [log(0.0565097)], [log(0.0565097)], [log(0.00722309)]]
    assert likelihood_parser(phyml_file_likelihoods) == likelihood_list

def test_error():
    with py.test.raises(IOError):
        likelihood_parser("no_file_name")
