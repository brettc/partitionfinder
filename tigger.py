import sys
from partfinder.alignment import Alignment
from partfinder._tiger import TigerDNA

if __name__ == "__main__":
    a = Alignment()
    a.read(sys.argv[1])
    tigger = TigerDNA()
    tigger.from_alignment(a)
    print tigger.calc_rates()

