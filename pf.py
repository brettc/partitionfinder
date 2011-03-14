import logging
log = logging.getLogger("main")

from optparse import OptionParser
import sys, os

from partfinder import config

def main():
    usage = """usage: python %prog [options] <foldername>

    PartitionFinder is designed to discover optimal partitioning schemes for
    DNA sequence alignments. It it also useful for finding the best model of
    sequence evolution for one or more partitions.

    The Input: <foldername>: the full path to a folder containing:
        - A configuration file (partition_finder.cfg)
        - A DNA alignment in Phylip format
    Take a look at the included 'example' folder for more details.

    The Output: A file in the same directory as the .cfg file, named
    'alignment_pf_output.txt' This file contains information on the best
    partitioning scheme, and the best model for each partiiton

    Usage Examples: 
        >python pf.py example
        Analyse what is in the 'example' sub-folder in the current folder.

        >python pf.py -v example
        Analyse what is in the 'example' sub-folder in the current folder, but
        show all the debug output

        >python pf.py -c ~/data/frogs
        Check the configuration files in the folder data/frogs in the current
        user's home folder.

        >python pf.py --force-restart ~/data/frogs
        Deletes any data produced by the previous runs (which is in
        ~/data/frogs/output) and starts afresh
    """
    parser = OptionParser(usage)
    parser.add_option(
        "-v", "--verbose",
        action="store_true", dest="verbose",
        help="show verbose (debug) output")
    parser.add_option(
        "-c", "--check-only",
        action="store_true", dest="check_only",
        help="just check the configuration files, don't do any processing")
    parser.add_option(
        "--force-restart",
        action="store_true", dest="force_restart",
        help="delete all previous output and start afresh (!)")
    
    options, args = parser.parse_args()

    # We should have one argument: the folder to read the configuration from
    if not args:
        # Otherwise exit, printing the help
        parser.print_help()
        return 2

    handler = logging.StreamHandler(sys.stdout)
    fmt = logging.Formatter('%(levelname)-8s | %(message)s')
    handler.setFormatter(fmt)
    logging.getLogger('').addHandler(handler)

    if options.verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO

    logging.getLogger().setLevel(logging.DEBUG)
    handler.setLevel(level)

    # Load, using the first argument as the folder
    try:
        config.initialise(args[0])
        config.load()
        # log.info("Configuration appears to be okay.")
        if options.check_only:
            log.info("Exiting without processing as requested...")
        else:
            # Now try processing everything....
            log.info("Beginning processing...")
            # config.process()
        # Successful exit
        log.info("Success: processing complete.")
        return 0
    except config.ConfigurationError:
        log.error("Configuration Failure: Please correct problems and rerun")
        # Any exceptions and we fail
    return 1


if __name__ == "__main__":
    # Well behaved unix programs exits with 0 on success...
    # sys.argv = ['arg', '-v', 'example']
    # sys.argv = ['arg', 'example']
    # sys.argv = ['arg', '-vc', 'example']
    sys.exit(main())


