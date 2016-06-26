# Notes on phyml and raxml

We have done our best to provide versions of PhyML and RAxML that will work for people on mac, windows, and linux. This document provides some guidance if the programs don't work for your particular setup.

##PhyML

We use the PhyML 3.1 binaries released here: http://www.atgc-montpellier.fr/download/binaries/phyml/PhyML-3.1.zip
We use the 64 bit linux binary.

If these don't work for you, try other binaries from that link, renaming them and replacing the version that comes with PF as instructed below.

If that doesn't work you will need to compile a version that works on your system. The source code is here: , the PhyML source code is here: https://github.com/stephaneguindon/phyml, use the PhyML docs for help, and use the PhyML google group if you get stuck.

Once you have found (or compiled) a version of PhyML that works on your system, double check that it works by running a test analysis, then rename it depending on your OS as follows:

mac: 'phyml'
windows: 'phyml.exe'
linux: 'phyml.linux'

then drop it into the /programs folder in partitionfinder, replacing the existing version.

## RAxML

The versions bundled in PF are from here:

https://github.com/stamatak/standard-RAxML/commit/3abe69bf5339ffe3e1d0abe9f6cc1e8980ab1be7

If the binaries that come with PF don't work for you, try other binaries (if you're on Windows) from the link above, renaming them and replacing the version that comes with PF as instructed below. Remember to get both the single-processor and PTHREADS versions.

### Building RAxML on a mac

cd to the directory where you downloaded the version of RAxML from the commit linked above

Then:

```
make -f Makefile.SSE3.PTHREADS.mac
rm *.o
make -f Makefile.SSE3.mac
rm *.o
```

Then rename the pthreads executable 'raxml_pthreads', and the non-pthreads executable 'raxml' and drop them in the programs folder. Overwite the old versions.

### Building RAxML on linux

For Linux, cd to the directory where you downloaded the version of RAxML from the commit linked above

```
make -f Makefile.SSE3.PTHREADS.gcc
rm *.o
make -f Makefile.SSE3.gcc
rm *.o
```

Then rename the pthreads executable 'raxml_pthreads.linux', and the non-pthreads executable 'raxml.linux' and drop them in the programs folder. Overwite the old versions.

### Building RAxML on Windows

This is hard, and if you want to do this you will need either expertise, or time and a lot of studying of the google group.

For our purposes, we used the pre-compiled executables that are in the version of RAxML from the commit linked above, which are in the folder within that repository called 'WindowsExecutables_v8.2.7_Alternate/', e.g. https://github.com/stamatak/standard-RAxML/tree/master/WindowsExecutables_v8.2.7_Alternate

If they don't work for you, try one of the other windows exectuble folders. Failing that, you will need to compile it yourself. Please don't request help with this on the PartitionFinder google group - we don't know how to do it.

Once you have a windows version to try, rename the pthreads executable 'raxml_pthreads.exe', and the non-pthreads executable 'raxml.exe' and drop them in the programs folder. Overwite the old versions.
