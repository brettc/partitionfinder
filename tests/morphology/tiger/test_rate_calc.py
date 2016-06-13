import pytest
import os
from partfinder import morph_tiger

def test_rates(set_parts):
    rates = morph_tiger.calculate_rates(set_parts)
    assert rates == [[0.375], [0.29166666666666663], [0.625], [0.5416666666666666], [1.0]]


if __name__ == '__main__':
    set_parts = [[[1,3],[2],[4]], [[1],[2],[3],[4]], [[1,2,3],[4]], [[1,2],[3,4]], [[1,2,3,4]]]
    test_rates(set_parts)
