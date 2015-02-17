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
*Partition Finder Morphology is in beta. Tests are passing, the software runs on test datasets. But it is still unreleased software. Please see the software terms of use for further info. See below for caveats on morphology.*
Currently, PartitionFinder works for morphology, when using the code in  branch feature/morphology2. To use this:

1. Clone the repository, as above.
2. Open Terminal (on a Mac) or Command Prompt (on Windows) and cd to the directory with PartitionFinder in it
3. Type 'git checkout feature/morphology2'
4. Now, PartitionFinder morphology can be executed. Similar to the nucelotide or protein versions, type:
```python
    python PartitionFinderMorphology.py morph/ cmdline-extras 'asc-corr=lewis'
```    
to execute the example, or replace <morph/> with your own folder containing a Phylip file and .cfg.

*Morphology has some special caveats.*
+ If you use automated partition discovery, this is based on the [TIGER](http://bioinf.nuim.ie/tiger/) of Cummins and McInerney (2011). Because this identifies characters that are dissimilar to other characters in the matrix (see paper for a discussion), some partitions may be quite small. That a partitioning scheme is supported statistically is not a guarantee that you will estimate a more correct topology using it, or that a Bayesian topology search will arrive at convergence using this scheme. We strongly suggest comparing trees between partitioned and unpartitioned runs, and, as always, _performing multiple rounds of topology estimation for any given set of parameters._
+ If all the data in your dataset are binary, specify 'binary' on line 10.
+ The command line option 'asc-corr=lewis' is an ascertainment bias correction. Unless you have collected invarient sites, you 
need to specify an ascertainment correction. 'asc-corr=lewis' is the correction described in Paul Lewis' [2001](http://sysbio.oxfordjournals.org/content/50/6/913) paper introducing the Mk model. More information can be found on 
the RAxML [website](http://sco.h-its.org/exelixis/resource/download/NewManual.pdf).
+ In order to take advantage of these important corrections, make sure you are using at least version 8.1.13
from the RAxML [github](https://github.com/stamatak/standard-RAxML/releases). 
+ Phylip is not a standard format for morphology, but it is very simple. It is simply the a one-line header with the name of species 
and characters, separated by a space. However, you do need to remove spaces from species names. If you would rather do this programmatically,
in the helper_scripts directory, find the script converter.py. It is called via:
```python
    python converter.py 'path to files' 'files extension' 'format of input files' 'format you need'
```
For example, to convert a nexus file in the current working directory, type:
```python
    python converter.py . .nex nexus phylip
```
This script depends on the Dendropy [library](https://pythonhosted.org/DendroPy/index.html).

For more details, read the manual.

# Special instructions for use of the iterative k-means algorithm

Follow the same instructions given in the manual, but with three important modifications:

1.  The new algorithm relies on some python dependencies. These can easily be downloaded and installed using Anaconda. You can [download Anaconda here](http://continuum.io/downloads).
2.  In the partitionfinder.cfg file, specify one data block for the entire alignment, e.g. for an alignment with 5,000 sites you would specify a single data block that reads, all = 1-5000;.
3.  In the partitionfinder.cfg file, modify the search to kmeans by specifying search = kmeans;

With these modifications, your analysis should be ready to go using the normal PartitionFinder commands.

Note that iterative kmeans has only been tested on Mac. It will also work on Linux, but the programs it is dependant on, i.e. raxml, phyml, and fast_TIGER need to be recompiled in Linux.