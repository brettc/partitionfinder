# PartitionFinder

[![Join the chat at https://gitter.im/brettc/partitionfinder](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/brettc/partitionfinder?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

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

* Make sure you have Python 2.7 installed first, if not, go to www.python.org/getit/

* For PartitionFinderProtein just substitute 'PartitionFinderProtein' for 'PartitionFinder' below

1.  Open Terminal (on a Mac) or Command Prompt (on Windows) and cd to the directory with PartitionFinder in it
2.  Run PartitionFinder by typing at the command prompt:

    python PartitionFinder.py example

This will run the included example analysis for PartitionFinder. More generally, the command line for PartitionFinder looks like this:

    python <PartitionFinder.py> <foldername>

where <PartitionFinder.py> is the full file-path to the PartitionFinder.py file
and <foldername> is the full filepath to a folder with a phylip alignemnt and associated .cfg file.

For more details, read the manual.
