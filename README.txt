# PartitionFinder
_________________

PartitionFinder is a Python program for simultaneously choosing partitioning
schemes and models of molecular evolution for DNA sequence data. It is useful
to use before running a phylogenetic analysis of DNA sequence data, in order
to decide how to divide up your sequence data into separate blocks before
analysis, and to simultaneously perform model selection on each of those
blocks.

# Operating System
PartitionFinder currently runs only on Mac and Windows.
All of the code was written with Linux in mind too, so if you are interested
in porting it to Linux, please get in touch (or just try it out!).

# Manual
is in the /docs folder. 

# Quick Start
# Make sure you have Python 2.7 installed first, if not, go to www.python.org/getit/

1.	Open Terminal (on a Mac) or Command Prompt (on Windows) and cd to the directory with PartitionFinder in it
2.	Run PartitionFinder by typing at the command prompt:

	python PartitionFinder.py example

    this will run the included example analysis for PartitionFinder
    
More generally, the command line for PartitionFinder looks like this:

pythong <PartitionFinder.py> <foldername>

where <PartitionFinder.py> is the full file-path to the PartitionFinder.py file
and <foldername> is the full filepath to a folder with a phylip alignemnt and associated .cfg file.

For more details, read the manual.