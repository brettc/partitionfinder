"""
Usage: PartitionFinder.py [options] <folder-path>

    PartitionFinder discovers optimal partitioning schemes and models for
    nucleotide and amino acid sequence alignments. See the examples below and
    the manual for more details.

Arguments:
    The <folder-path> must contain:
        - A configuration file named "partition_finder.cfg"
        - A nucleotide/aa alignment in Phylip format
Output:
    All output files can be found in the "analysis" folder, which can be
    found under the supplied <folder-name>.

Options:
  -c, --check-only      Just check the configuration files, don't do any processing
  -f, --force-restart   Delete all previous output and start afresh (!)
  -h, --help            Show this help message and exit
  -p, --processes <N>   Number of concurrent processes to use [default: all]
  -q, --quick           Avoid anything slow (like writing schemes at each step), useful for very large datasets
  -r, --raxml           Use RAxML, rather than PhyML, to do the analysis
  -v, --verbose         Show all debug logging information

Advanced Options:
  --weights <O:B:M:A>   Assign different weights to (O)verall rate for a
                        subset, the (B)ase/amino acid frequencies, (M)odel
                        parameters, and (A)lpha value [default: 1:0:0:0].
  --kmeans-opt <N>      Choose which kmeans algorithm to use [default: 1]:
                            1: Use site likelihoods only (PhyML and RAxML)
                            2: Use site rates only (PhyML only)
                            3: Use site likelihoods and site rates (PhyML only)
                            4: Use site likelihoods from gamma rate categories
                               (PhyML only)
  --rcluster-percent <N>
                        Defines the percentage of schemes that relaxed
                        clustering algorithm will consider before it stops
                        looking [default: 10].
  --rcluster-max <N>
                        Defines the number of schemes that relaxed
                        clustering algorithm will consider before it stops
                        looking [default: all].
  --cmdline-extras <command-string>
                        Add additional commands to the phyml or raxml command
                        (may need to be quoted)

Debugging Options:
  --debug-output <region:region:...>
                        Show debug information for only specified code regions
  --list-regions        List possible debug regions (implies --check-only)
  --show-python-exceptions
                        If errors occur, print the python exceptions
  --save-phylofiles     Save all of the phyml or raxml output. This takes a
                        lot of space(!)
  --dump-results        Dump final results to a binary file.
  --compare-results     Compare the results to previously dumped binary file.

Examples:
    > python PartitionFinder.py example
    Analyse what is in the 'example' sub-folder in the current folder.

    > python PartitionFinder.py --force-restart ~/data/frogs
    Deletes any data produced by the previous runs (which is in
    ~/data/frogs/output) and starts afresh
"""

# TODO add some of these advanced examples
import os
import threadpool

advance = """
Advanced Examples:

    lines that PartitionFinder uses.  e.g. if you want to change the accuracy of lnL
    calculations ('-e' option in raxml), or use
    multi-threaded versions of raxml that require you to
    specify the number of threads you will let raxml use ('-T' option in raxml.  E.g. you might specify this: --cmndline_extras ' -e 2.0 -T 10 ' N.B. MAKE SURE YOU PUT YOUR EXTRAS IN
        QUOTES, and only use this command if you really know what you're doing and are very familiar with raxml
        and PartitionFinder

                        This will affect how
                        subsets are clustered together. For instance:
                            --cluster_weights '1, 2, 5, 1', would weight the
                            base freqeuncies 2x more than the overall rate, the
                            model parameters 5x more, and the alpha parameter
                            the same as the model rate
"""
import docopt
import version

class Options(object):
    """All program options get stored here, along with their validation

    Some come from the command line. Others from the config.
    """

    def __init__(self):
        # NOTE: we let most defaults get populated by docopt processing
        pass

    def parse_args(self, args):
        opts = docopt.docopt(__doc__, version=version.get_version(),
                                  argv=args)
        self.validate(opts)

    def validate(self, opts):
        """Process the dictionary returned by docopt"""
        for option, value in opts.items():
            if option.startswith('--'):
                option = option[2:]
            elif option.startswith("<"):
                option = option[1:-1]
            option = option.replace('-', '_')

            checker = getattr(self, "check_{}".format(option), None)
            if checker:
                # The checker can transform the value too
                value = checker(value)

            if hasattr(self, option):
                # This can only happen if someone has screwed up the
                # configuration. It means that there are two options with the
                # same name
                raise KeyError

            setattr(self, option, value)

    # Define checkers for options here. Just make a function called
    # check_OPTION_NAME, and return the value it should hold, or raise an error
    def check_folder_path(self, f):
        return f
        # if not os.path.exists(f):
        #     raise ""

    def check_processes(self, v):
        if v == 'all':
            v = threadpool.get_cpu_count()
        else:
            v = int(v)
            if v < 1: v = 1
        return v

