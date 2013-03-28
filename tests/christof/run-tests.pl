#!/usr/bin/perl

### In order to add new test search for __DATA__ toward the end of this file.


use strict;

my $res;
my $count_pass = 0;
my $count_fail = 0;
my $count_expected_fail = 0;
my $count_new_reference = 0;
my $count_new_reference_failed = 0;
my $count_dies = 0;

my @commands = <DATA>;
my $com;
foreach $com (@commands)
{
    ## Remove leading and trailing spaces:
    chomp($com);
    $com =~ s/^\s+//;
    $com =~ s/\s+$//;

    if (!($com eq "") && !($com =~ /^#/))
    {
	print "$com\n";

	$res = system($com);
    
##	print "$res\n";
   
	if ($res == 0)
	{
	    ++$count_pass;
	}
	elsif ($res == 256) ## exit(1)
	{
	    ++$count_expected_fail;
	}
	elsif ($res == 512) ## die or exit(2)
	{
	    ++$count_dies;
	}
	elsif ($res == 768) ## exit(3)
	{
	    ++$count_new_reference_failed;
	}
	elsif ($res == 1024) ## exit(4)
	{
	    ++$count_new_reference;
	}
	else ## exit(-1)
	{
	    ++$count_fail;
	}
    }
}

print "Number of new reference analysis created:                                      $count_new_reference\n";
print "Number of new reference analysis failed:                                       $count_new_reference_failed\n";
print "Number of tests passed:                                                        $count_pass\n";
print "Number of tests made to fail and they did fail with the correct error message: $count_expected_fail\n"; 
print "Number of tests failed unexpectedly:                                           $count_fail\n";
print "Number of tests in which test skript failed (bad parameters etc.)              $count_dies\n";



## Add new tests here:
## After the __DATA__ keywords, each line defines a test.
## These lines are passed to the single-pf-test.pl skript.
###################
## Entries must be:
###################
## - Script name: single-pf-test.pl
## - Path to the Partitionfinder.py or PartitionFinderProtein.py file
## - Folder with cfg and alignment file
## - A string that is expected in the log file. E.g. "Processing complete" or "The 'search = clustering' option is only availalbe when using raxml"
#    In this way we can run tests which have to fail and we check for the correct error message to occur.
## - Parameters to PartitionFinder, e.g. "--raxml"      




__DATA__
########### Basic tests using example data sets (6)
#single-pf-test.pl /Daten/PartitionFinder-Versions/pf/partitionfinder-clustering3_v2/PartitionFinderProtein.py /Volumes/Hugo/PF-tests/test-aaExample-greedy-raxml "Processing complete" "--raxml"
#single-pf-test.pl /Daten/PartitionFinder-Versions/pf/partitionfinder-clustering3_v2/PartitionFinderProtein.py /Volumes/Hugo/PF-tests/test-aaExample-greedy-phyml "Processing complete"
#single-pf-test.pl /Daten/PartitionFinder-Versions/pf/partitionfinder-clustering3_v2/PartitionFinder.py /Volumes/Hugo/PF-tests/test-nucExample-greedy-raxml       "Processing complete" "--raxml"
#single-pf-test.pl /Daten/PartitionFinder-Versions/pf/partitionfinder-clustering3_v2/PartitionFinder.py /Volumes/Hugo/PF-tests/test-nucExample-greedy-phyml       "Processing complete"
#single-pf-test.pl /Daten/PartitionFinder-Versions/pf/partitionfinder-clustering3_v2/PartitionFinder.py /Volumes/Hugo/PF-tests/test-nucExample-clustering-phyml   "The 'search = clustering' option is only availalbe when using raxml"
#single-pf-test.pl /Daten/PartitionFinder-Versions/pf/partitionfinder-clustering3_v2/PartitionFinder.py /Volumes/Hugo/PF-tests/test-nucExample-clustering-raxml   "Processing complete" "--raxml"

########### Testing cfg file parsing: (15)
#single-pf-test.pl /Daten/PartitionFinder-Versions/pf/partitionfinder-clustering3_v2/PartitionFinderProtein.py /Volumes/Hugo/PF-tests/test-aaExample-schemeBrackets-greedy-phyml   "It looks like the '[schemes]' option might be missing or in the wrong place"
#single-pf-test.pl /Daten/PartitionFinder-Versions/pf/partitionfinder-clustering3_v2/PartitionFinderProtein.py /Volumes/Hugo/PF-tests/test-aaExample-filespec-in-2-lines-greedy-phyml   "Processing complete"
#single-pf-test.pl /Daten/PartitionFinder-Versions/pf/partitionfinder-clustering3_v2/PartitionFinderProtein.py /Volumes/Hugo/PF-tests/test-aaExample-umlinked-greedy-phyml   "The only valid options for 'branchlengths' are: 'linked', 'unlinked'"
#single-pf-test.pl /Daten/PartitionFinder-Versions/pf/partitionfinder-clustering3_v2/PartitionFinderProtein.py /Volumes/Hugo/PF-tests/test-aaExample-order-model-selection-greedy-phyml " It looks like the 'branchlengths' option might be missing or in the wrong place"
#single-pf-test.pl /Daten/PartitionFinder-Versions/pf/partitionfinder-clustering3_v2/PartitionFinderProtein.py /Volumes/Hugo/PF-tests/test-aaExample-model-selection-greedy-phyml "The line causing the problem is this: 'model_selection = AIC, BIC;'"
#single-pf-test.pl /Daten/PartitionFinder-Versions/pf/partitionfinder-clustering3_v2/PartitionFinderProtein.py /Volumes/Hugo/PF-tests/test-aaExample-data-blocks1-greedy-phyml "Expected \"=\" (at char"
#single-pf-test.pl /Daten/PartitionFinder-Versions/pf/partitionfinder-clustering3_v2/PartitionFinderProtein.py /Volumes/Hugo/PF-tests/test-aaExample-data-blocks2-greedy-phyml "It looks like the '[schemes]' option might be missing or in the wrong place"

#single-pf-test.pl /Daten/PartitionFinder-Versions/pf/partitionfinder-clustering3_v2/PartitionFinderProtein.py /Volumes/Hugo/PF-tests/test-aaExample-data-blocks3-greedy-phyml "It looks like the '[schemes]' option might be missing or in the wrong place"
#single-pf-test.pl /Daten/PartitionFinder-Versions/pf/partitionfinder-clustering3_v2/PartitionFinderProtein.py /Volumes/Hugo/PF-tests/test-aaExample-data-blocks4-greedy-phyml "Partition(EF1a, 625-949\1, 1-1\1, 3-3\1) overlaps with previously defined partitions"
#single-pf-test.pl /Daten/PartitionFinder-Versions/pf/partitionfinder-clustering3_v2/PartitionFinderProtein.py /Volumes/Hugo/PF-tests/test-aaExample-data-blocks5-greedy-phyml "Partition 'EF1a' has internal overlap"
#single-pf-test.pl /Daten/PartitionFinder-Versions/pf/partitionfinder-clustering3_v2/PartitionFinderProtein.py /Volumes/Hugo/PF-tests/test-aaExample-data-blocks6-greedy-phyml "Expected W:(0123...)"
### First data blocks starts at 0. No error is raised. But the created data sets are wrong. The first position is doubled. The total number of AAs is increased by one.
#single-pf-test.pl /Daten/PartitionFinder-Versions/pf/partitionfinder-clustering3_v2/PartitionFinderProtein.py /Volumes/Hugo/PF-tests/test-aaExample-data-blocks7-greedy-phyml "Processing complete"
#single-pf-test.pl /Daten/PartitionFinder-Versions/pf/partitionfinder-clustering3_v2/PartitionFinderProtein.py /Volumes/Hugo/PF-tests/test-aaExample-data-blocks8-greedy-phyml "Site 1000 is specified in [data_blocks], but the alignment only has 949 sites. Please check."
### Order of data-blocks is irrelevant
#single-pf-test.pl /Daten/PartitionFinder-Versions/pf/partitionfinder-clustering3_v2/PartitionFinderProtein.py /Volumes/Hugo/PF-tests/test-aaExample-data-blocks9-greedy-phyml "Processing complete"
### Non optimal error message
#single-pf-test.pl /Daten/PartitionFinder-Versions/pf/partitionfinder-clustering3_v2/PartitionFinderProtein.py /Volumes/Hugo/PF-tests/test-aaExample-data-blocks10-greedy-phyml "It looks like the '[schemes]' option might be missing or in the wrong place"



########### Testing model specification/parsing: (8)
#single-pf-test.pl /Daten/PartitionFinder-Versions/pf/partitionfinder-clustering3_v2/PartitionFinderProtein.py /Volumes/Hugo/PF-tests/test-aaExample-modellistAll-greedy-phyml  "Processing complete" 
#single-pf-test.pl /Daten/PartitionFinder-Versions/pf/partitionfinder-clustering3_v2/PartitionFinderProtein.py /Volumes/Hugo/PF-tests/test-aaExample-modellist-error1-greedy-phyml  "'WAG+I+G+I' is not a valid model for phylogeny program phyml." 
#single-pf-test.pl /Daten/PartitionFinder-Versions/pf/partitionfinder-clustering3_v2/PartitionFinderProtein.py /Volumes/Hugo/PF-tests/test-aaExample-modellist-error2-greedy-phyml  "Processing complete"
#single-pf-test.pl /Daten/PartitionFinder-Versions/pf/partitionfinder-clustering3_v2/PartitionFinderProtein.py /Volumes/Hugo/PF-tests/test-aaExample-modellist-error3-greedy-phyml  "'WAG+LG+F+I' is not a valid model for phylogeny program phyml." 
#single-pf-test.pl /Daten/PartitionFinder-Versions/pf/partitionfinder-clustering3_v2/PartitionFinderProtein.py /Volumes/Hugo/PF-tests/test-aaExample-modellist-error4-greedy-raxml  "'WAG+F' is not a valid model for phylogeny program raxml."  "--raxml"
#single-pf-test.pl /Daten/PartitionFinder-Versions/pf/partitionfinder-clustering3_v2/PartitionFinderProtein.py /Volumes/Hugo/PF-tests/test-aaExample-modellist-error5-greedy-raxml  "The line causing the problem is this: 'models = MtArt+I+G+F WAG+I+G;'"  "--raxml"
#single-pf-test.pl /Daten/PartitionFinder-Versions/pf/partitionfinder-clustering3_v2/PartitionFinderProtein.py /Volumes/Hugo/PF-tests/test-aaExample-modellist-error6-greedy-raxml   "'MtArt+I+G+F' is not a valid model for phylogeny program raxml"  "--raxml"
### Testing multiline models specification:
#single-pf-test.pl /Daten/PartitionFinder-Versions/pf/partitionfinder-clustering3_v2/PartitionFinderProtein.py /Volumes/Hugo/PF-tests/test-aaExample-modellist-error7-greedy-phyml   "Processing complete"

########## Testing sequence length: (6)
#single-pf-test.pl /Daten/PartitionFinder-Versions/pf/partitionfinder-clustering3_v2/PartitionFinderProtein.py /Volumes/Hugo/PF-tests/test-aaExample-seqLen99-greedy-raxml "Processing complete" "--raxml"
#single-pf-test.pl /Daten/PartitionFinder-Versions/pf/partitionfinder-clustering3_v2/PartitionFinderProtein.py /Volumes/Hugo/PF-tests/test-aaExample-seqLen99-greedy-phyml "Processing complete"
### Maximum length of sequences that are distinuished in raxml is 99. If two sequences differ at posiiton 100 they are not distinguished.
#single-pf-test.pl /Daten/PartitionFinder-Versions/pf/partitionfinder-clustering3_v2/PartitionFinderProtein.py /Volumes/Hugo/PF-tests/test-aaExample-seqLen100-greedy-raxml "Sequence names of taxon 1 and 2 are identical" "--raxml"
#single-pf-test.pl /Daten/PartitionFinder-Versions/pf/partitionfinder-clustering3_v2/PartitionFinderProtein.py /Volumes/Hugo/PF-tests/test-aaExample-seqLen100-greedy-phyml "Processing complete"
### Misleading error message 
#single-pf-test.pl /Daten/PartitionFinder-Versions/pf/partitionfinder-clustering3_v2/PartitionFinderProtein.py /Volumes/Hugo/PF-tests/test-aaExample-seqLen101-greedy-raxml "A common cause of this error is having whitespace," "--raxml"
### Misleading error message
#single-pf-test.pl /Daten/PartitionFinder-Versions/pf/partitionfinder-clustering3_v2/PartitionFinderProtein.py /Volumes/Hugo/PF-tests/test-aaExample-seqLen101-greedy-phyml "A common cause of this error is having whitespace,"


