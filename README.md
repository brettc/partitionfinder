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

#Morphology notes.

Currently, PartitionFinder works for morphology, when using the code in  branch feature/morphology2. To use this:

1. Clone the repository, as above.
2. Open Terminal (on a Mac) or Command Prompt (on Windows) and cd to the directory with PartitionFinder in it
3. Type 'git checkout feature/morphology2'
4. Now, PartitionFinder morphology can be executed. Similar to the nucelotide or protein versions, type:

    python <PartitionFinderMorphology.py> <morph/> <--cmdline-extras 'asc-corr=lewis'>
    
to execute the example, or replace <morph/> with your own folder containing a Phylip file and .cfg.

*Morphology has some special caveats.*
+ If all the data in your dataset are binary, specify 'binary' on line 10.
+ The command line option 'asc-corr=lewis' is an ascertainment bias correction. Unless you have collected invarient sites, you 
need to specify an ascertainment correction. 'asc-corr=lewis' is the correction described in Paul Lewis' [2001](http://sysbio.oxfordjournals.org/content/50/6/913) paper introducing the Mk model. More information can be found on 
the RAxML [website](http://sco.h-its.org/exelixis/resource/download/NewManual.pdf).
+ In order to take advantage of these important corrections, make sure you are using at least version 8.1.13
from the RAxML [github](https://github.com/stamatak/standard-RAxML/releases). 
+ Phylip is not a standard format for morphology, but it is very simple. It is simply the a one-line header with the name of species 
and characters, separated by a space. However, you do need to remove spaces from species names. If you would rather do this programmatically,
in the helper_scripts directory, find the script converter.py. It is called via:

    python <converter.py> <'path to files' 'extension of files' 'format of input files' 'format you'd like exported'>

This script depends on the Dendropy [library](https://pythonhosted.org/DendroPy/index.html).

For more details, read the manual.
