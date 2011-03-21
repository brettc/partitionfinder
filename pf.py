import logging, os
curdir, f = os.path.split(__file__)
config_path = os.path.join(curdir, 'logging.cfg')
from logging import config as _logconfig
_logconfig.fileConfig(config_path)

log = logging.getLogger("main")

from optparse import OptionParser
import sys

from partfinder import config, analysis

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
    parser.add_option(
        "-p", "--processes", 
        type="int", dest="processes", default=1, metavar="N",
        help="Number of concurrent processes to use."
        " Use -1 to match the number of cpus on the machine."
        " The default is to use 1."
    )
    
    options, args = parser.parse_args()

    # We should have one argument: the folder to read the configuration from
    if not args:
        # Otherwise exit, printing the help
        parser.print_help()
        return 2

    # handler = logging.StreamHandler(sys.stdout)
    # fmt = logging.Formatter('%(levelname)-8s | %(message)s')
    # handler.setFormatter(fmt)
    # logging.getLogger('').addHandler(handler)

    # if options.verbose:
        # level = logging.DEBUG
    # else:
        # level = logging.INFO

    # logging.getLogger().setLevel(logging.DEBUG)
    # handler.setLevel(level)

    # Load, using the first argument as the folder
    try:
        config.load(args[0])
        if options.check_only:
            log.info("Exiting without processing (because of the -c/--check-only option ...")
        else:
            # Now try processing everything....
            s = config.settings
            anal = analysis.Analysis(
                s.alignment_path,
                s.analysis_path,
                options.force_restart,
                threads=options.processes,
            )
            if s.search_algorithm == 'all':
                anal.analyse_all_possible(s.models)
            elif s.search_algorithm == 'user':
                anal.analyse_current_schemes(s.models)
            else:
                log.error("Search algorithm %s is not yet implemented", 
                          s.search_algorithm)
                raise NotImplemented

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
    sys.exit(main())


