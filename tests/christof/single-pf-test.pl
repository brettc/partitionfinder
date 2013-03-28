#!/usr/bin/perl

## They should all be standard modules according to http://de.selfhtml.org/perl/module/standardmodule.htm
use strict;
use POSIX qw(strftime);
use File::Copy;
use File::Compare;
use File::Path 'rmtree';

########## Global variable that need to be set ##############


my $pf_args_default = "-v --save-phylofiles --debug-output=raxml,analysis,subset";
my $expected_fail = 0;


# Linux and mac
# On Unix machines we supress most of the output.
my $todevnull = "";

if (-e "/dev/null")
{
    $todevnull = " 1> /dev/null 2> /dev/null";
}


sub last_line_of_file
{
    my $filename = $_[0];

    my $pos=-1;
    my $char;
    my $already_nonblank = 0;

    open IN, "$filename" or die "Could not open file $filename which is required for this test.";

    while (seek (IN,$pos--,2))
    {
	read IN,$char,1;
	last if ($char eq "\n" and $already_nonblank == 1);
	$already_nonblank = 1 if ($char ne "\n");
    }
    
    my $last_line = <IN>;
    close IN;
    chomp($last_line);
    return $last_line;
}

sub string_occurrs_in_file
{
    my $filename = $_[0];
    my $str      = $_[1];

    open IN, "$filename" or die "Could not open file $filename which is required for this test.";

    $str = quotemeta($str);

#    $str =~ s/\[//g;
#    $str =~ s/\]//g;
#    $str =~ s/\+//g;
#    print "$str\n";

    while (<IN>) ## for all lines in this file
    {
#	s/\[//g;
#	s/\]//g;
#	s/\+//g;
#	print "$_\n";

	my $line = $_;

#	print "$line\n";

	if ($line =~ m/$str/) ## Check this line matches $str
	{
	    return 1;
	}
    }
    return 0; ## Nothing found
}


##
## We only call this function before conducting a test.
## If we clean up afterwords, we erase all information we need to find out what went wrong.
##
sub clean_up_analysis
{
    my $data = $_[0];

    if (-e "$data/log.txt")
    {
	print "Removing log.txt file from previous run.\n";
	unlink("$data/log.txt");
    }
    if (-d "$data/analysis")
    {
	print "Removing analysis directory from previous run.\n";
	rmtree([ "$data/analysis" ]);
    }
}




# Windows:
# my $devnull = "NUL";


my $argc = @ARGV;

if ($argc <3 || $argc > 4 )
{
    die "Error: single-pf-test.pl requires three or four parameters.\nUsage: single-pf-test.pl partition-finder-version directory-with-data-and-cfg-file \"Expected log file string such as Processing complete\" [partition-finder-arguments].\n"
}

my $pf                       = $ARGV[0];
my $data                     = $ARGV[1];
my $expected_log_file_string = $ARGV[2];

my $pf_args = "";

if ($argc > 3)
{
   $pf_args = $ARGV[3];
}

if (! -e $pf)
{
    die "TEST SKRIPT PARAMETER ERROR: The program \"$pf\" specified as first argument to this script does not exists.\n";
}

if (! -d $data)
{
    die "TEST SKRIPT PARAMETER ERROR: The data direcotry \"$data\" specified as second argument to this script does not exists.\n";
}

my $isref = 0;

if (! (-d "$data/analysis_ref") ) # We have no reference yet, so we create one:
{
    $isref = 1;
}


##########################
### Start analysis:
##########################

## Before we start we clean up any previous data:
clean_up_analysis("$data");

my $pf_command = "python $pf $pf_args_default $pf_args $data $todevnull";

print "Starting partition finder with the following command line:\n";
print "$pf_command\n";
print "This may take a while ...\n";

system($pf_command);



####################################################
### Did we succeed? Do some basic tests:
####################################################

if (! (-e "$data/log.txt") )
{
    print "TEST FAILED: No log file has been created.\n\n";
    exit(-1);
}

my $expected_string_found_in_new_logfile  = string_occurrs_in_file("$data/log.txt", "$expected_log_file_string");

if ($expected_string_found_in_new_logfile == 0 && ! -d "$data/analysis")
{
    print "TEST FAILED: Unexpectedly, no analysis folder has been created and the log file does not contain the expected string.\n";
    print "Details can be found in the log file: $data/log.txt\n\n"; 
    exit(-1);
}

if ($expected_string_found_in_new_logfile == 0 && $isref ) 
{
    print "SETUP OF REFENENCE ANALYSIS FAILED: Expected string not found in log file.\n";
    print "Details can be found in the log file: $data/log.txt\n\n";
    exit(3);
}



## Handle important special case: For some errors we will encounter, no analysis foulder will be created.
## In this case we create the "analysis" folder to store the log.txt file.
## The "analysis" folder will be moved to "analysis_ref" later.

if ($expected_string_found_in_new_logfile == 1 &&  !(-d "$data/analysis") )
{
    mkdir("$data/analysis");

    if (!$isref)  ## If this is not a refence analysis, it will be an analysis that is expected to fail
    {
	$expected_fail = 1;
    }
}

###########################
### Move log and rename result folder:
###########################

my $timestamp = strftime "%Y.%m.%d_%H.%M.%S", localtime;

# Remove trailing "/" if it exists:
$data =~ s/\/$//;

##my $com1 = "mv $data/log.txt ${data}analysis";
##my $com2 = "mv $data/analysis $data/analysis_$timestamp";

move("$data/log.txt", "${data}/analysis") or die "Failed to move ${data}/log.txt\n";

if ($isref) # We have no reference yet, so we create one:
{
    print("Analysis $data/analysis will be used as new reference anaylsis.\n");
    move("$data/analysis", "$data/analysis_ref") or die "Failed to move $data/analysis\n\n";
    exit(4);
}
else
{
    move("$data/analysis", "$data/analysis_$timestamp") or die "Failed to move $data/analysis\n";
}

###########################
### Check result:
###########################

if ($expected_string_found_in_new_logfile == 0)
{
    if (!$isref)
    {
	print "TEST FAILED: Expected string not found in $data/analysis_$timestamp/log.txt\n\n";
	exit(-1);
    }
    else
    {
	print "FAILED TO CREATE REFERENCE ANALYSIS: Expected string not found in $data/analysis_ref/log.txt\n\n";
	exit(3);
    }
}

## From this point on, we can be sure that ($expected_string_found_in_new_logfile == 1)!!!


if (! $isref)
{
    my $expected_string_found_in_ref_logfile  = string_occurrs_in_file("$data/analysis_ref/log.txt", $expected_log_file_string);

    if ($expected_string_found_in_ref_logfile == 0)
    {
	print("TEST SETUP PROBLEM: Expected string found in new but not in refence log file.\nCheck: $data/analysis_ref/log.txt\nand\n$data/analysis_$timestamp/log.txt.\n\n");
	exit(3);
    }
}

## In case we did not find the expected string in either the new or reference log file we already exited.
## So far everything looks good. But let us do two final tests:


if (!(-e "$data/analysis_ref/best_schemes.txt") || !("$data/analysis_ref/all_schemes.txt") )
{
    ## We can only get to this point if the expected string has been found in the log file.
    ## Now if the best_schemes.txt or the all_schemes.txt file in the reference analysis are missing
    ## this has be an expected fail
    $expected_fail = 1;
}


if ($expected_fail == 1)
{
    print "Test passed. PartitionFinder had to tail for $data.\n\n";
    exit(1);
}



if (compare("$data/analysis_ref/best_schemes.txt","$data/analysis_$timestamp/best_schemes.txt") != 0)
{
    print ("TEST FAILED: Unequal best_schemes.txt files for test $data\n");
    print "Please compare the log files: $data/analysis_ref/log.txt\nand\n$data/analysis_$timestamp/log.txt.\n\n"; 
    exit(-1);
}

if (compare("$data/analysis_ref/all_schemes.txt","$data/analysis_$timestamp/all_schemes.txt") != 0)
{
    print ("TEST FAILED: Unequal all_schemes.txt files for test $data\n");
    print "Please compare the log files: $data/analysis_ref/log.txt\nand\n$data/analysis_$timestamp/log.txt.\n\n"; 
    exit(-1);
}

print "Test passed for $data.\n\n";

exit(0);

## In principle one could go into more details if the run has to fail and to check for the occurrence of the correct error message.
