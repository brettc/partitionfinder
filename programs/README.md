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

The versions bundled in PF are from here: https://github.com/stamatak/standard-RAxML/commit/222e44ae7133ee1039b85853ffc5fb6964ae9c9e (* if you're not used to GitHub, this points to the specific version of the code that we used). There are lots of varieties of RAxML you can use. But for PartitionFinder you need the single-processor version, and (if you want to use the --ml-tree in PF) the PTHREADS version. Specifically, the versions we use are:

mac: Makefile.SSE3.mac; Makefile.SSE3.PTHREADS.mac
linux: Makefile.SSE3.gcc; Makefile.SSE3.PTHREADS.gcc (compiled on Ubuntu)
Windows: /WindowsExecutables_v8.2.4_Alternate/raxmlHPC-SSE3.exe; /WindowsExecutables_v8.2.4_Alternate/raxmlHPC-PTHREADS-SSE3.exe

If the binaries that come with PF don't work for you, try other binaries (if you're on Windows) from the link above, renaming them and replacing the version that comes with PF as instructed below. Remember to get both the single-processor and PTHREADS versions.

If the other binaries don't work, or you're not using Windows, you will need to compile a version that works on your system. Follow the link to the source code above (you could also try the lastest version, but we can't guarantee that the latest version will work with PartitionFinder), use the RAxML docs to help you compile both the single processor and PTHREADS binaries, use the RAxML google group if you get stuck.

Once you have found (or compiled) the two versions of RAxML you need, double check that both work by running test analyses, then rename them depending on your OS as follows:

Single processor
mac: 'raxml'
windows: 'raxml.exe'
linux: 'raxml.linux'

PTHREADS
mac: 'raxml_pthreads'
windows: 'raxml_pthreads.exe'
linux: 'raxml_pthreads.linux'

then drop both of these into the /programs folder in partitionfinder, replacing the existing versions.