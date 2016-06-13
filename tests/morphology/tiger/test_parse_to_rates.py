import pytest
import os
from partfinder import morph_tiger
from partfinder.alignment import Alignment

MORPH_DATA = """
5 25
tax1        1101100000111111021200001
tax2        110110010?0011100212?0000
tax3        ?11210?1010?00101000?010?
tax4        1021011201000010111111111
tax5        10?1011311000120121211111

"""

MORPH_ALIGN = Alignment()
MORPH_ALIGN.parse(MORPH_DATA)

def test_set_parts():
    set_parts = morph_tiger.create_set_parts(MORPH_ALIGN)
    assert set_parts == [[[0, 1, 3, 4]], [[3, 4], [0, 1, 2]], [[0, 1], [2], [3]], [[0, 1, 3, 4], [2]], [[3, 4], [0, 1, 2]], [[0, 1, 2], [3, 4]], [[0, 1], [3, 4]], [[0], [1, 2], [3], [4]], [[0, 1, 2, 3], [4]], [[0], [2, 3, 4]], [[1, 2, 3, 4], [0]], [[1, 3, 4], [0]], [[2, 3, 4], [0, 1]], [[2, 3], [0, 1, 4]], [[0, 1, 2, 3], [4]], [[1, 2, 3, 4], [0]], [[0, 1], [2, 3, 4]], [[2], [3], [0, 1, 4]], [[2], [0, 1, 3, 4]], [[2], [3], [0, 1, 4]], [[0], [3, 4]], [[0, 1, 2], [3, 4]], [[0, 1], [2, 3, 4]], [[0, 1, 2], [3, 4]], [[1], [0, 3, 4]]]

def test_axpi():
    set_parts = morph_tiger.create_set_parts(MORPH_ALIGN)
    axpi_0 = morph_tiger.axpi(set_parts[0], set_parts[1])
    assert axpi_0 == 0.5
    axpi_1 = morph_tiger.axpi(set_parts[0], set_parts[2])
    two_thirds = 2.0/3
    assert axpi_1 == two_thirds

def test_rates():
    set_parts = morph_tiger.create_set_parts(MORPH_ALIGN)
    rates = morph_tiger.calculate_rates(set_parts)
    print rates

test_rates()
