from analysis import Analysis
from partition import all_partitions

def show_partitions():
    for p in all_partitions:
        print p

if __name__ == "__main__":
    show_partitions()