# PartitionFinder

PartitionFinder and PartitionFinderProtein are Python programs for simultaneously 
choosing partitioning schemes and models of molecular evolution for sequence data. 
You can use them before running a phylogenetic analysis, in order
to decide how to divide up your sequence data into separate blocks before
analysis, and to simultaneously perform model selection on each of those
blocks.

# Operating System

Mac and Windows are supported.
All of the code was written with Linux in mind too, so if you are interested
in porting it to Linux, please get in touch (or just try it out!).

# Manual

is in the /docs folder. 

# Quick Start

* Make sure you have the Anaconda Python distribution version 2.3.0 or higher, and make sure it's the one that ships with Python 2, and that the version number of Python is at least 2.7.10 (you do not want the anaconda distribution that gives you Python 3.x). http://continuum.io/downloads.

* For PartitionFinderProtein just substitute 'PartitionFinderProtein' for 'PartitionFinder' below

1.  Open Terminal (on a Mac) or Command Prompt (on Windows) and cd to the directory with PartitionFinder in it
2.  Run PartitionFinder by typing at the command prompt:

    python PartitionFinder.py examples/nucleotide

This will run the included example analysis for PartitionFinder. More generally, the command line for PartitionFinder looks like this:

    python <PartitionFinder.py> <foldername>

where <PartitionFinder.py> is the full file-path to the PartitionFinder.py file
and <foldername> is the full filepath to a folder with a phylip alignemnt and associated .cfg file.

For more details, read the manual.
