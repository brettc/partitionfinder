import pytest
from partfinder import alignment


def test_mat():
    '''This test should pass'''
    test = """
10 2
Allosaurus_fragilis                 11
Sinraptor                           11
Dilong_paradoxus                    11
Eotyrannus_lengi                    11
Tyrannosaurus_rex                   11
Gorgosaurus_libratus                11
Tanycolagreus_topwilsoni            11
Coelurus_fragilis                   11
Ornitholestes_hermanni              10
Huaxiagnathus_orientalis            10

    """
    alignment.AlignmentParser(test)		

