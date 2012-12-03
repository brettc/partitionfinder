#Copyright (C) 2012 Robert Lanfear and Brett Calcott
#
#This program is free software: you can redistribute it and/or modify it
#under the terms of the GNU General Public License as published by the
#Free Software Foundation, either version 3 of the License, or (at your
#option) any later version.
#
#This program is distributed in the hope that it will be useful, but
#WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#General Public License for more details. You should have received a copy
#of the GNU General Public License along with this program.  If not, see
#<http://www.gnu.org/licenses/>. PartitionFinder also includes the PhyML
#program, the RAxML program, the PyParsing library, and the python-cluster library
#all of which are protected by their own licenses and conditions, using
#PartitionFinder implies that you agree with those licences and conditions as well.

import logging
import sys

logging.basicConfig(
    format="%(levelname)-8s | %(asctime)s | %(message)s",
    level=logging.INFO
)

# curdir = os.path.dirname(os.path.abspath(__file__))
# rootdir, here = os.path.split(curdir)
# config_path = os.path.join(rootdir, 'logging.cfg')
# from logging import config as _logconfig
# _logconfig.fileConfig(config_path)

log = logging.getLogger("main")
from optparse import OptionParser

# We import everything here as it forces all of debug regions to be loaded
import version
import config
import analysis_method
import util
import reporter
import progress
import datetime
import parser
import raxml
import phyml


def debug_arg_callback(option, opt, value, parser):
    setattr(parser.values, option.dest, value.split(','))


def get_debug_regions():
    mlogger = logging.Logger.manager
    return mlogger.loggerDict.keys()


def set_debug_regions(regions):
    if regions is None:
        return
    valid_regions = set(get_debug_regions())
    if 'all' in regions:
        regions = valid_regions
    else:
        regions = set(regions)
        errors = set()
        for r in regions:
            if r not in valid_regions:
                log.error("'%s' is not a valid debug region", r)
                errors.add(r)
        if errors:
            raise util.PartitionFinderError

    for r in regions:
        logging.getLogger(r).setLevel(logging.DEBUG)

    # Enhance the format
    fmt = logging.Formatter("%(levelname)-8s | %(asctime)s | %(name)-10s | %(message)s")
    logging.getLogger("").handlers[0].setFormatter(fmt)


def main(name, datatype):
    v = version.get_git_version()
    log.info("------------- %s %s -----------------", name, v)
    start_time = datetime.datetime.now().replace(microsecond=0)  # start the clock ticking
    usage = """usage: python %prog [options] <foldername>

    PartitionFinder and PartitionFinderProtein are designed to discover optimal
    partitioning schemes for nucleotide and amino acid sequence alignments.
    They are also useful for finding the best model of sequence evolution for datasets.

    The Input: <foldername>: the full path to a folder containing:
        - A configuration file (partition_finder.cfg)
        - A nucleotide/aa alignment in Phylip format
    Take a look at the included 'example' folder for more details.

    The Output: A file in the same directory as the .cfg file, named
    'analysis' This file contains information on the best
    partitioning scheme, and the best model for each partiiton

    Usage Examples:
        >python %prog example
        Analyse what is in the 'example' sub-folder in the current folder.

        >python %prog -v example
        Analyse what is in the 'example' sub-folder in the current folder, but
        show all the debug output

        >python %prog -c ~/data/frogs
        Check the configuration files in the folder data/frogs in the current
        user's home folder.

        >python %prog --force-restart ~/data/frogs
        Deletes any data produced by the previous runs (which is in
        ~/data/frogs/output) and starts afresh
    """
    parser = OptionParser(usage)
    parser.add_option(
        "-v", "--verbose",
        action="store_true", dest="verbose",
        help="show debug logging information (equivalent to --debug-out=all)")
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
        type="int", dest="processes", default=-1, metavar="N",
        help="Number of concurrent processes to use."
        " Use -1 to match the number of cpus on the machine."
        " The default is to use -1.")
    parser.add_option(
        "--show-python-exceptions",
        action="store_true", dest="show_python_exceptions",
        help="If errors occur, print the python exceptions")
    parser.add_option(
        "--save-phylofiles",
        action="store_true", dest="save_phylofiles",
        help="save all of the phyml or raxml output. This can take a lot of space(!)")
    parser.add_option(
        "--dump-results",
        action="store_true", dest="dump_results",
        help="Dump all results to a binary file. "
        "This is only of use for testing purposes.")
    parser.add_option(
        "--compare-results",
        action="store_true", dest="compare_results",
        help="Compare the results to previously dumped binary results. "
        "This is only of use for testing purposes.")
    parser.add_option(
        "--raxml",
        action="store_true", dest="raxml",
        help="Use RAxML (rather than PhyML) to do the analysis. See the manual"
    )
    parser.add_option(
        "--cmdline-extras",
        type="str", dest="cmdline_extras", default="", metavar="N",
        help="Add additional commands to the phyml or raxml commandlines that PF uses."
        "This can be useful e.g. if you want to change the accuracy of lnL calculations"
        " ('-e' option in raxml), or use multi-threaded versions of raxml that require"
        " you to specify the number of threads you will let raxml use ('-T' option in "
        "raxml. E.g. you might specify this: --cmndline_extras ' -e 2.0 -T 10 '"
        " N.B. MAKE SURE YOU PUT YOUR EXTRAS IN QUOTES, and only use this command if you"
        " really know what you're doing and are very familiar with raxml and"
        " PartitionFinder"
    )
    parser.add_option(
        "--cluster-weights",
        type="str", dest="cluster_weights", default=None, metavar="N",
        help="Mainly for algorithm development. Only use it if you know what you're doing."
        "A list of weights to use in the clustering algorithm. This list allows you "
        "to assign different weights to the overall rate for a subset, the base/amino acid "
        "frequencies, and the model parameters. This will affect how subsets are "
        "clustered together. For instance: --cluster_weights '1, 2, 5', would weight "
        "the base freqeuncies 2x more than the overall rate, and the model parameters 5x "
        "more."
    )
    parser.add_option(
        '--debug-output',
        type='string',
        action='callback',
        dest='debug_output',
        metavar="REGION,REGION,...",
        callback=debug_arg_callback,
        help="(advanced option) Provide a list of debug regions to output extra "
        "information about what the program is doing."
        " Possible regions are 'all' or any of {%s}."
        % ",".join(get_debug_regions())
    )

    options, args = parser.parse_args()

    #default to phyml
    if options.raxml == 1:
        options.phylogeny_program = 'raxml'
    else:
        options.phylogeny_program = 'phyml'

    # Error checking
    if options.dump_results and options.compare_results:
        log.error("You can't dump and compare results in one run!")
        log.error("Please select just one of these")

    # We should have one argument: the folder to read the configuration from
    if not args:
        # Otherwise exit, printing the help
        parser.print_help()
        return 2

    #before we start, let's check the python version is above 2.7 but lower than 3.0
    python_version = float(
        "%d.%d" % (sys.version_info[0], sys.version_info[1]))

    log.info("You have Python version %.1f" % python_version)

    if python_version < 2.7:
        log.error("Your Python version is %.1f, but this program requires Python 2.7. "
                  "Please upgrade to version 2.7 by visiting www.python.org/getit, or by following"
                  " the instructions in the PartitionFinder manual." % python_version)
        return 0

    if python_version > 3.0:
        log.warning("Your Python version is %.1f. This program was not built to run with "
                    "version 3 or higher. To guarantee success, please use Python 2.7.x" % python_version)

    # Load, using the first argument as the folder
    try:
        if options.verbose:
            set_debug_regions(['all'])
        else:
            set_debug_regions(options.debug_output)
        cfg = config.Configuration(datatype, options.phylogeny_program,
                                   options.save_phylofiles, options.cmdline_extras, options.cluster_weights)
        # Set up the progress callback
        p = progress.TextProgress(cfg)
        cfg.load_base_path(args[0])

        if options.check_only:
            log.info("Exiting without processing (because of the -c/--check-only option ...")
        else:


            # Now try processing everything....
            method = analysis_method.choose_method(cfg.search)
            reporter.TextReporter(cfg)
            anal = method(cfg,
                          options.force_restart,
                          options.processes)
            results = anal.analyse()

            if options.dump_results:
                results.dump(cfg)
            elif options.compare_results:
                results.compare(cfg)

        # Successful exit
        end_time = datetime.datetime.now().replace(microsecond=0)
        processing_time = end_time - start_time

        log.info("Total processing time: %s (h:m:s)" % processing_time)
        log.info("Processing complete.")

        return 0

    except util.PartitionFinderError:
        log.error("Failed to run. See previous errors.")
        if options.show_python_exceptions:
            raise

    except KeyboardInterrupt:
        log.error("User interrupted the Program")

    return 1
