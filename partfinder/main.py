# Copyright (C) 2012 Robert Lanfear and Brett Calcott
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details. You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# PartitionFinder also includes the PhyML program, the RAxML program, and the
# PyParsing library, all of which are protected by their own licenses and
# conditions, using PartitionFinder implies that you agree with those licences
# and conditions as well.

__VERSION__ = "2.0.0"

import logging
import sys
import shlex
import logtools

logging.basicConfig(
    format="%(levelname)-8s | %(asctime)s | %(message)s",
    level=logging.INFO
)

log = logtools.get_logger()
from optparse import OptionParser

# We import everything here as it forces all of debug regions to be loaded
import config
import analysis_method
import util
import reporter
import progress
import datetime
import parser
import raxml
import phyml


def debug_arg_callback(option, opt, value, theparser):
    setattr(theparser.values, option.dest, value.split(','))


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
            return errors

    for r in regions:
        logging.getLogger(r).setLevel(logging.DEBUG)

    # Enhance the format
    fmt = logging.Formatter(
        "%(levelname)-8s | %(asctime)s | %(name)-10s | %(message)s")
    logging.getLogger("").handlers[0].setFormatter(fmt)

    return None

def parse_args(datatype, cmdargs=None):
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
    op = OptionParser(usage)
    op.add_option(
        "-v", "--verbose",
        action="store_true", dest="verbose",
        help="show debug logging information (equivalent to --debug-out=all)")
    op.add_option(
        "-c", "--check-only",
        action="store_true", dest="check_only",
        help="just check the configuration files, don't do any processing")
    op.add_option(
        "-f", "--force-restart",
        action="store_true", dest="force_restart",
        help="delete all previous output and start afresh (!)")
    op.add_option(
        "-p", "--processes",
        type="int", dest="processes", default=-1, metavar="N",
        help="Number of concurrent processes to use."
             " Use -1 to match the number of cpus on the machine."
             " The default is to use -1.")
    op.add_option(
        "--show-python-exceptions",
        action="store_true", dest="show_python_exceptions",
        help="If errors occur, print the python exceptions")
    op.add_option(
        "--save-phylofiles",
        action="store_true", dest="save_phylofiles",
        help="save all of the phyml or raxml output. This can take a lot of space(!)")
    op.add_option(
        "--dump-results",
        action="store_true", dest="dump_results",
        help="Dump all results to a binary file. "
             "This is only of use for testing purposes.")
    op.add_option(
        "--compare-results",
        action="store_true", dest="compare_results",
        help="Compare the results to previously dumped binary results. "
             "This is only of use for testing purposes.")
    op.add_option(
        "-q", "--quick",
        action="store_true", dest="quick", default=False,
        help="Avoid anything slow (like writing schemes at each step)," 
             "useful for very large datasets."
    )
    op.add_option(
        "-r", "--raxml",
        action="store_true", dest="raxml",
        help="Use RAxML (rather than PhyML) to do the analysis. See the manual"
    )

    op.add_option(
        "-n", "--no-ml-tree",
        action="store_true", dest="no_ml_tree", default=False,
        help="Estimate a starting tree with NJ (PhyML) or MP (RaxML) instead"
             " of the default which is to estimate a starting tree with ML "
             " using in RAxML. Not recommended."
    )

    op.add_option(
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
    op.add_option(
        "--weights",
        type="str", dest="cluster_weights", default=None, metavar="N",
        help="Mainly for algorithm development. Only use it if you know what you're doing."
             "A list of weights to use in the clustering algorithms. This list allows you "
             "to assign different weights to: the overall rate for a subset, the base/amino acid "
             "frequencies, model parameters, and alpha value. This will affect how subsets are "
             "clustered together. For instance: --cluster_weights '1, 2, 5, 1', would weight "
             "the base freqeuncies 2x more than the overall rate, the model parameters 5x "
             "more, and the alpha parameter the same as the model rate"
    )
    op.add_option(
        "--kmeans",
        type="str", dest="kmeans", default='entropy', metavar="type",
        help="This defines which sitewise values to use: entropy or tiger "
             "\n--kmeans entropy: use entropies for sitewise values"
             "\n--kmeans tiger: use TIGER rates for sitewise values (only valid for Morphology)"
    )
    op.add_option(
        "--rcluster-percent",
        type="float", dest="cluster_percent", default=10.0, metavar="N",
        help="This defines the proportion of possible schemes that the relaxed clustering"
             " algorithm will consider before it stops looking. The default is 10%."
             "\ne.g. --rcluster-percent 10.0"
    )
    op.add_option(
        "--rcluster-max",
        type="int", dest="cluster_max", default=-987654321, metavar="N",
        help="This defines the number of possible schemes that the relaxed clustering"
             " algorithm will consider before it stops looking. The default is to look at "
             "the larger value out of 1000, and 10 times the number of data blocks you have."
             "\ne.g. --rcluster-max 1000"
    )
    op.add_option(
        "--min-subset-size",
        type="int", dest="min_subset_size", default=False, metavar="N",
        help="This defines the minimum subset size that the kmeans and rcluster"
             " algorithm will accept. Subsets smaller than this "
             " will be merged at with other subsets at the end of the algorithm"
             " (for kmeans) or at the start of the algorithm (for rcluster)."
             " See manual for details. The default value for kmeans is 100."
             " The default value for rcluster is to ignore this option."
             "\ne.g. --min-subset-size 100"
    )
    op.add_option(
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
    op.add_option(
        "--all-states",
        action="store_true", dest="all_states", default=False,
        help="In the kmeans and rcluster algorithms, this stipulates that PartitionFinder "
             "should not produce subsets that do not have all possible states present. E.g."
             " for DNA sequence data, all subsets in the final scheme must have A, C, T, "
             " and G nucleotides present. This can occasionally be useful for downstream "
             " analyses, particularly concerning amino acid datasets."
    )

    op.add_option(
        '--profile',
        action="store_true",
        help="Output profiling information after running (this will slow everything down!)")

    if cmdargs is None:
        options, args = op.parse_args()
    else:
        options, args = op.parse_args(cmdargs)

    options.datatype = datatype
    # We should have one argument: the folder to read the configuration from
    if not args:
        op.print_help()
    else:
        check_options(op, options)

    return options, args


def check_options(op, options):
    # Error checking
    if options.dump_results and options.compare_results:
        op.error(
            "options --dump_results and --compare_results are mutually exclusive!")

    if options.verbose:
        set_debug_regions(['all'])
    else:
        errors = set_debug_regions(options.debug_output)
        if errors is not None:
            bad = ",".join(list(errors))
            op.error("Invalid debug regions: %s" % bad)

    # Default to phyml
    if options.raxml == 1:
        options.phylogeny_program = 'raxml'
    else:
        options.phylogeny_program = 'phyml'

    #A warning for people using the Pthreads version of RAxML
    if options.cmdline_extras.count("-T") > 0:
        log.warning("""
            It looks like you're using a Pthreads version of RAxML. Be aware
            that the default behaviour of PartitionFinder is to run one version
            of RAxML per available processor. This might not be what you want
            with Pthreads, since the minimum number of threads per RAxML run is
            2 (i.e. -T 2). Make sure to limit the total number of RAxML runs
            you start using the -p option in PartitionFinder.  Specifically,
            the total number of processors you will use with the Pthreads
            version is the number you set via the -T option in
            --cmdline-extras, multiplied by the number of processors you set
            via the -p option in PartitionFinder.  You should also be aware
            that the Pthreads version of RAxML has a rare but known bug on some
            platforms. This bug results in infinite likelihood values if it
            happens on your dataset, PartitionFinder will give an error. In
            that case you should switch back to using a single-threaded version
            of RAxML, e.g. the SSE3 or AVX version.  See the manual for more
            info.""")


def check_python_version():
    """Check the python version is above 2.7 but lower than 3.0"""

    python_version = float(
        "%d.%d" % (sys.version_info[0], sys.version_info[1]))

    log.info("You have Python version %.1f" % python_version)

    if python_version < 2.7:
        log.error("""
            Your Python version is %.1f, but this program requires Python 2.7.
            Please upgrade to version 2.7 by visiting www.python.org/getit,
            or  by following the instructions in the PartitionFinder manual.
            """ % python_version)
        return 0

    if python_version > 3.0:
        log.warning("""
            Your Python version is %.1f. This program was not built to run
            with version 3 or higher. To guarantee success, please use
            Python 2.7.x""" % python_version)

def run_analysis(cfg, options):
    # Now try processing everything....
    method = analysis_method.choose_method(cfg.search)
    reporter.TextReporter(cfg)
    anal = method(cfg, options.force_restart, options.processes)
    results = anal.analyse()
    if options.dump_results:
        results.dump(cfg)
    elif options.compare_results:
        results.compare(cfg)

def profile_analysis(cfg, options):
    import cProfile, pstats
    cProfile.runctx('run_analysis(cfg, options)', globals(), locals(), filename='profile.output')
    p = pstats.Stats('profile.output')
    p.sort_stats('time').print_stats(20)
    p.sort_stats('cumtime').print_stats(20)
    # p.strip_dirs().sort_stats(-1).print_stats()

def main(name, datatype, passed_args=None):

    # If passed_args is None, this will use sys.argv
    options, args = parse_args(datatype, passed_args)
    if not args:
        # Help has already been printed
        return 2

    log.info("------------- %s %s -----------------" % (name, __VERSION__))
    start_time = datetime.datetime.now().replace(
        microsecond=0)  # start the clock ticking

    check_python_version()

    if passed_args is None:
        cmdline = " ".join(sys.argv)
    else:
        cmdline = " ".join(passed_args)

    log.info("Command-line arguments used: %s" % cmdline)

    # Load, using the first argument as the folder
    try:
        # TODO: just pass the options in!
        config.the_config.init(datatype,
                                   options.phylogeny_program,
                                   options.save_phylofiles,
                                   options.cmdline_extras,
                                   options.cluster_weights,
                                   options.cluster_percent,
                                   options.cluster_max,
                                   options.kmeans, 
                                   options.quick,
                                   options.min_subset_size,
                                   options.all_states,
                                   options.no_ml_tree)
        cfg = config.the_config

        # Set up the progress callback
        progress.TextProgress(cfg)
        cfg.load_base_path(args[0])

        if options.check_only:
            log.info("""
            Exiting without processing (because of the
            -c/--check-only option ...
            """)
        else:
            try:
                if options.profile:
                    profile_analysis(cfg, options)
                else:
                    run_analysis(cfg, options)
            finally:
                # Make sure that we reset the configuration
                cfg.reset()

        # Successful exit
        end_time = datetime.datetime.now().replace(microsecond=0)
        processing_time = end_time - start_time

        log.info("Total processing time: %s (h:m:s)" % processing_time)
        log.info("Processing complete.")

        return 0

    except util.ExternalProgramError as e:
        log.error("""A program that Partitionfinder uses failed. Output
                  follows, in case it's helpful for finding the problem""")
        log.error("%s", e.stdout)
        log.error("%s", e.stderr)
        if options.show_python_exceptions or passed_args is not None:
            raise

    except util.ParseError:
        log.error("""We failed to parse the output of an external program (See
                  previous errors).  This is probably due to an update to the
                  program which has changed the output.""")
        if options.show_python_exceptions or passed_args is not None:
            raise

    except util.PartitionFinderError:
        log.error("Failed to run. See previous errors.")
        # Reraise if we were called by call_main, or if the options is set
        if options.show_python_exceptions or passed_args is not None:
            raise

    except (KeyboardInterrupt, SystemExit):
        log.error("User interrupted the Program")

    return 1


def call_main(datatype, cmdline):
    cmdargs = shlex.split(cmdline)
    main("<main.py:call_main(%s)>" % datatype, datatype, cmdargs)
