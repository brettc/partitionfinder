import sys
from partfinder.alignment import Alignment
from partfinder._tiger import TigerDNA
from pathlib import Path
import numpy

if __name__ == "__main__":
    a = Alignment()
    filepath = Path(sys.argv[1])
    a.read(str(filepath))

    tigger = TigerDNA()
    tigger.build_bitsets(a)
    rates = tigger.calc_rates()
    output_path = str(filepath.with_suffix('.tigger8'))
    numpy.savetxt(output_path, rates, fmt="%5f", delimiter='\n')


    

