#Using PartitionFinder for enormous datasets

##Install Anaconda
This just makes sure that you have the latest Python libraries that PartitionFinder needs.

1. Go to this address and follow the instructions for your system:
http://docs.continuum.io/anaconda/install.html

2. Check that it worked by typing this at the command prompt:
    conda

It should print out a lot of help stuff for anaconda.

##Install PartitionFinder
We need the development version of PartitionFinder (it will be PF2, but it's not quite ready for general release, yet):

1. Go to here to get the zip file:
https://github.com/brettc/partitionfinder/archive/develop.zip

2. Unzip develop.zip which will create a folder called partitionfinder-develop

3. If you are using Linux, you will need to download and compile your own version of RAxML.
3.1 Download the PF version of RAxML from here:
https://github.com/brettc/standard-RAxML/archive/master.zip

3.2 Compile using the Makefile.SSE3.gcc, or the Makefile.AVX.gcc version, depending on your chip architecture, following the instructions here (do not use the PTHREADs or MPI versions, and do not use the RAxML from Alexis Statmatakis' github account):
https://github.com/brettc/standard-RAxML

3.3 Rename the compiled executable "raxml"

3.4 In the 'partitionfinder-develop/programs' folder, delete 'raxml'
copy-paste your newly compiled 'raxml' into partitionfinder-develop/programs folder


##Run PartitionFinder
Setting up the .cfg file still works exactly as described in the manaul (in the /docs folder). There are just a few additional options you should know about, that are not yet described in the manual.


--rcluster-max
defines the number of partitioning scheme that PF will explore at each step of the analysis. More is better, and the default is 1000. As a rough rule of thumb, I suggest using 2x the number of data blocks in your .cfg file. If your run is too slow, make --rcluster-max smaller. We have to be pragmatic, after all. But if your run finishes quickly, try running it again after doubling --rcluster max. 

-q
this stops PF from writing unnecessary files. And helps you to stay friends with the folks that manage the clusters you might use. 

By default PF will use ALL of your CPUs. To control the number of cpus, use -p:

-p 10
This will use 10 cpus (or less than 10 if you have less than 10).

##An example
Here's a worked example of analysing a massive dataset with PartitionFinder.

###Example partition_finder.cfg file 
Let's imagine we're running a very large amino acid dataset, with 1000's of species, and 1000's of data blocks.  Our .cfg file might look something like this:

    alignment = alignment.phy;

    branchlengths = linked;

    models = LG+G, LG+G+F;

    model_selection = AICc;

    [data_blocks]

    Gene1 = 1 - 609 ;
    Gene2 = 610 - 1409;
    Gene3 = 1410 - 1560;
    ... [LOTS MORE DATA BLOCKS HERE]

    [schemes]
    search = rcluster;

Let's go through each line in turn. 

'alignment' just points the alignment file. Name it whatever your alignment is called.

'branchlengths' should be set to linked. You can try 'unlinked' if you like, but this invovles estimating a LOT more paratmeters and is almost never preferred over linked branch lengths. In almost all cases, it's fine to just assume that 'linked' is the best option.

'models' defines the models we want to use. See the manual for a list of what's possible (we have to stick to raxml models here). I've defined a very small model list here, because analysing very big datasets can take a long time. So it's worth thinking carefully about how to reduce the number of models you're analysing. In general, twice the number of models will take twice as long to analyse. Also be aware that, if you're planning to analyse your data in RAxML after you do the partitioning, your model list should not mix models with and without +I (read the RAxML manual to see why - RAxML can't mix these kinds of model in a single analysis).

'model_selection': just choose your favourite metric from AICc and BIC. The AIC is not recommended, becuase it's asymptotically the same as the AICc as dat blocks get very big, but the AICc is more sensible for smaller data blocks.

'[data blocks]': look at the manual to see how to define these.

'search': stick to 'rcluster' for giant datasets. Just use --rcluster-max at the command line to control the algorithm (see above).

###Example command line
I'm assuming you've read the manual before you read this. I'm not going to explain anything that's already in the manual.

Here's an example of running PF the simplest way we can for our .cfg file above (assuming there are no errors in your .cfg file or alignment):

    python ~/blah/PartitionFinderProtein.py ~/blah/analysis --raxml

We need the "--raxml" to tell PF to use raxml - the rcluster algorithm only works with raxml. And in general if you're working with massive datasets, you'll be using raxml too.

It's probably a good idea with massive datasets to avoid any unneccesary file writing by using the -q option:

    python ~/blah/PartitionFinderProtein.py ~/blah/analysis --raxml -q

If we want to control the relaxed clustering algorithm, we have a few options. We can make it more thorough than the default by increasing the number of schemes it looks at:

    python ~/blah/PartitionFinderProtein.py ~/blah/analysis --raxml -q --rcluster-max 5000

If you want to control the number of processors, use -p:

    python ~/blah/PartitionFinderProtein.py ~/blah/analysis --raxml -q --rcluster-max 5000 -p 20

And if you want to control the way that similarity is defined among subsets, you can use the "--weights" option (see the manual for how this works). For example, we could assign equal weight to rates, amino acid frequencies, and the gamma parameter, but ignore model parameters (sensible for an amino acid dataset because they are not estimated from the data):

    python ~/blah/PartitionFinderProtein.py ~/blah/analysis --raxml -q --rcluster-max 5000 -p 20 --weights "1, 1, 0, 1"



  --weights <O:B:M:A>   Assign different weights to (O)verall rate for a
                        subset, the (B)ase/amino acid frequencies, (M)odel
                        parameters, and (A)lpha value [default: 1:0:0:0].








