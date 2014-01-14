#Using PartitionFinder for enormous datasets

##Install Anaconda
This just makes sure that you have the latest Python libraries that PartitionFinder needs.

1. Go to this address and follow the instructions for your system:
http://docs.continuum.io/anaconda/install.html

2. Check that it worked by typing this at the command prompt:
    > conda

It should print out a lot of help stuff for anaconda.

##Install PartitionFinder
We need the development version of PartitionFinder (it will be PF2, but it's not quite
ready for general release, yet):

Go to here to get the zip file:
https://github.com/brettc/partitionfinder/archive/develop.zip

Then unzip develop.zip which will create a folder called partitionfinder-develop

If you are using Linux, you will need to download and compile your own version of RAxML.
First, download the PF version of RAxML from here:
https://github.com/brettc/standard-RAxML/archive/master.zip

Compile using the Makefile.SSE3.gcc, or the Makefile.AVX.gcc version, depending on your 
chip architecture, following the instructions here (do not use the PTHREADs or MPI 
versions, and do not use the RAxML from Alexis Statmatakis' github account):
https://github.com/brettc/standard-RAxML

Rename the compiled executable "raxml"

In the 'partitionfinder-develop/programs' folder, delete 'raxml'
copy-paste your newly compiled 'raxml' into partitionfinder-develop/programs folder


PF Run Setup
************
Make your .cfg file with the following settings (you need to define your own data blocks!):

branchlengths = linked;
models = LG+G, LG+G+F, LG4M+G, LG4M+G+F, JTT+G, JTT+G+F;
model_selection = AICc;
search = rcluster;

Read the manual in partitionfinder-develop/docs if you're not sure what to do. 

Run PartitionFinderProtein as described in the manual. Here's an example:

python ~blah/PartitionFinderProtein.py ~/blah/partition_finder.cfg --raxml --rcluster-max 1000 -q

That commandline contains two useful options that are not listed in the manual:

--rcluster-max
defines the number of partitioning scheme that PF will explore at each step of the analysis
more is better. The default is 1000, but a rough rule of thumb is to use 2x the number of 
data blocks in your .cfg file. If your run is too slow, make --rcluster-max smaller. We 
have to be pragmatic, after all.

-q
this stops PF from writing unnecessary files. And helps you to stay friends with the folks
that manage the clusters you might use. 

By default PF will use ALL of your CPUs. To control the number of cpus, use -p:

-p 10
This will use 10 cpus (or less than 10 if you have less than 10).



