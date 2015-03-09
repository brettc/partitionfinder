import sys
from partfinder.alignment import Alignment
from partfinder._tiger import TigerDNA
from pathlib import Path
import numpy
import pandas as pd

if __name__ == "__main__":
    a = Alignment()
    filepath = Path(sys.argv[1])
    a.read(str(filepath))

    tigger = TigerDNA()
    tigger.build_bitsets(a)
    rates = tigger.calc_rates()
    output_path = str(filepath.with_suffix('.tigger'))
    numpy.savetxt(output_path, rates, fmt="%5f", delimiter='\n')
    x = pd.DataFrame(tigger.calc_array())
    output_path = str(filepath.with_suffix('.tigger_array'))
    x.to_pickle(str(output_path))


    

