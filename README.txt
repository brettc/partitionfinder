# PartitionFinder
_________________

PartitionFinder is a Python program for simultaneously choosing partitioning
schemes and models of molecular evolution for DNA sequence data. It is useful
to use before running a phylogenetic analysis of DNA sequence data, in order
to decide how to divide up your sequence data into separate blocks before
analysis, and to simultaneously perform model selection on each of those
blocks.

# Operating System
PartitionFinder currently runs only on Intel Macs. Sorry.
All of the code was written with Windows in mind too, so if you are interested
in porting it to Windows, it shouldn't be too much work.

# Manual
is in the /docs folder. 

# Quick Start

1.	Open Terminal and cd to the directory with PartitionFinder in it
2.	Run PartitionFinder by typing at the command prompt:

	python PartitionFinder.py <foldername>

	where ‘<foldername>’ is the full file path to a folder containing a 
	Phylip alignment and a PartitionFinder configuration file called 
	‘partition_finder.cfg’. E.g.

	python PartitionFinder.py /Users/Rob/new_analysis