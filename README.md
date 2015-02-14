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

# Special instructions for use of the iterative k-means algorithm

Follow the same instructions given in the manual, but with two important modifications:

1.  In the partitionfinder.cfg file, specify one data block for the entire alignment, e.g. for an alignment with 5,000 sites you would specify a single data block that reads, all = 1-5000;.
2.  In the partitionfinder.cfg file, modify the search to kmeans by specifying search = kmeans;

With these modifications, your analysis should be ready to go using the normal PartitionFinder commands.

Note that iterative kmeans has only been tested on Mac. It will also work on Linux, but the programs it is dependant on, i.e. raxml, phyml, and fast_TIGER need to be recompiled in Linux.