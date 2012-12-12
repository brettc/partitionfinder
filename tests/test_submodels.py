from partfinder.submodels import get_submodels, count_all_schemes

def test_consistency():
    known_results = [
        [0, 0, 0, 0],
        [0, 0, 0, 1],
        [0, 0, 1, 0],
        [0, 0, 1, 1],
        [0, 0, 1, 2],
        [0, 1, 0, 0],
        [0, 1, 0, 1],
        [0, 1, 0, 2],
        [0, 1, 1, 0],
        [0, 1, 1, 1],
        [0, 1, 1, 2],
        [0, 1, 2, 0],
        [0, 1, 2, 1],
        [0, 1, 2, 2],
        [0, 1, 2, 3],
    ]
    submodels = get_submodels(4)
    assert submodels ==  known_results

def test_scheme_lengths():
    assert count_all_schemes(1) == 1
    assert count_all_schemes(5) == 52
    assert count_all_schemes(10) == 115975
