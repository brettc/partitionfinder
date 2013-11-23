import logging
log = logging.getLogger("config")

import os
import fnmatch
import scheme
import subset
import parser
import util
import progress

class ConfigurationError(util.PartitionFinderError):
    pass


class Configuration(object):
    # List of valid options. The first one is the default
    options = {
        'branchlengths': ['linked', 'unlinked'],
        'model_selection': ['aic', 'aicc', 'bic'],
        'search': ['all', 'user', 'greedy', 'hcluster', 'rcluster', 'kmeans',
            'kmeans_wss', 'kmeans_greedy']
    }

    def __init__(self, datatype="DNA", phylogeny_program='phyml',
                 save_phylofiles=False, cmdline_extras="", cluster_weights=None,
                 cluster_percent=10, kmeans_opt=1):

        log.info("------------- Configuring Parameters -------------")
        # Only required if user adds them
        self.user_schemes = scheme.SchemeSet()
        self.user_subsets = []
        self.user_subsets_by_name = {}

        self.save_phylofiles = save_phylofiles
        self.progress = progress.NoProgress(self)
        self.cmdline_extras = cmdline_extras
        self.cluster_percent = float(cluster_percent)
        self.kmeans_opt = kmeans_opt

        # Record this
        self.base_path = '.'
        self.alignment = None
        self.user_tree = None
        self.old_working_directory = None

        # Some basic checking of the setup, so that we don't hit too many
        # problems later
        if datatype != "DNA" and datatype != "protein":
            log.error("datatype must be 'DNA' or 'protein'")
            raise ConfigurationError

        log.info("Setting datatype to '%s'", datatype)
        self.datatype = datatype

        if phylogeny_program != "phyml" and phylogeny_program != "raxml":
            log.error("Phylogeny program must be 'phyml' or 'raxml'")
            raise ConfigurationError

        # Import the right processor
        self.processor = __import__(phylogeny_program.lower(), globals())

        log.info("Setting phylogeny program to '%s'", phylogeny_program)
        self.phylogeny_program = phylogeny_program

        # Now find the phylogeny_program before we do any changing of dirs
        self.find_programs()

        if cluster_weights is None:
            # default weights - just use overall rates of subsets. Based on
            # 2013 analyses.
            self.cluster_weights = {"rate": 1, "freqs": 0,
                                    "model": 0, "alpha": 0}
        else:
            # TODO. Is there a more robust way to do this...
            # Brett say "YES. But this will do for now..."
            cluster_weights = [x.strip() for x in cluster_weights.split(",")]

            # now we check that it's a list of exactly four numbers
            if len(cluster_weights) != 4:
                log.error("Your --cluster_weights argument should have exactly 4"
                          " numbers separated by commas, but it has %d ('%s') "
                          "Please check and try again", len(cluster_weights), cluster_weights)
                raise ConfigurationError

            for thing in cluster_weights:
                try:
                    num = float(eval(thing))
                    assert num >= 0
                except:
                    log.error("Unable to understand your --cluster_weights argument."
                              " It should look like this: --cluster_weights '1,2,3,6'. "
                              "Please double check that you included quotes, "
                              "and four numbers greater than or equal to zero "
                              "separated by commas. Then try again. "
                              "The part that I couldn't understand is this: '%s'" % thing)
                    raise ConfigurationError

            log.info("Setting cluster_weights to: "
                     "subset_rate = %s, freqs = %s, model = %s, alpha %s"
                     % (cluster_weights[0], cluster_weights[1],
                        cluster_weights[2], cluster_weights[3]))

            self.cluster_weights = {}
            self.cluster_weights["rate"] = float(eval(cluster_weights[0]))
            self.cluster_weights["freqs"] = float(eval(cluster_weights[1]))
            self.cluster_weights["model"] = float(eval(cluster_weights[2]))
            self.cluster_weights["alpha"] = float(eval(cluster_weights[3]))

        # Set the defaults into the class. These can be reset by calling
        # set_option(...)
        for o, v in self.options.items():
            # Could call self.set_option here -- but it might confuse users
            setattr(self, o, v[0])

        try:
            assert self.cluster_percent >= 0.0
            assert self.cluster_percent <= 100.0
        except:

            log.error("The rcluster-percent variable must be between 0.0 to 100.0, yours "
                      "is %.2f. Please check and try again." % self.cluster_percent)
            raise ConfigurationError

        log.debug("Setting rcluster-percent to %.2f" % self.cluster_percent)

        if kmeans_opt < 1 or kmeans_opt > 4:
            log.error(
                "The --kmeans-opt setting must be 1, 2, 3, or 4. Please check and restart")
            raise ConfigurationError

    def find_programs(self):
        pth = os.path.abspath(__file__)
        # Split off the name and the directory...
        pth, notused = os.path.split(pth)
        pth, notused = os.path.split(pth)
        pth = os.path.join(pth, "programs")
        pth = os.path.normpath(pth)
        self.program_path = pth

        # TODO This is bullshit---Need to make the config global
        util.program_path = pth
        log.info("Program path is here %s", self.program_path)

    def reset(self):
        if self.old_working_directory is not None:
            log.debug(
                "Returning to original path: %s", self.old_working_directory)
            os.chdir(self.old_working_directory)
        log.debug(
            "Cleaning out all subsets (There are %d)...", subset.count_subsets())
        subset.clear_subsets()

    def find_config_file(self, pth):
        """Try and get the base folder and config file from the path"""
        if os.path.isfile(pth):
            # Is it a config file
            pth, ext = os.path.splitext(pth)
            folder, filename = os.path.split(pth)
            filename += ext
            if ext == '.cfg':
                return folder, filename
            # We still need a filename
        else:
            folder = pth

        util.check_folder_exists(folder)

        # Now let's find the filename. Just return the first hit.
        for filename in os.listdir(folder):
            if fnmatch.fnmatch(filename, '*.cfg'):
                return folder, filename

        log.error("Cannot find a configuration file in "
                  "working folder '%s'", folder)

        raise ConfigurationError

    def load_base_path(self, pth):
        """Load using a base path folder"""
        # Allow for user and environment variables
        pth = os.path.expanduser(pth)
        pth = os.path.expandvars(pth)
        pth = os.path.normpath(pth)

        folder, filename = self.find_config_file(pth)

        # check that user didn't enter a file instead of a folder
        # if os.path.isfile(pth):
            # log.error("The second argument of the commandline currently
            # points to a file, but it should point to the folder that contains
            # the alignment and .cfg files, please check.")
            # raise ConfigurationError

        self.set_base_path(folder)

        # From now on we refer to relative paths
        config_path = os.path.join(self.base_path, filename)
        log.debug("About to search for partition_finder.cfg file...")
        config_path = os.path.join(self.base_path, "partition_finder.cfg")
        util.check_file_exists(config_path)

        self._output_folders = []
        self.register_output_folders()

        self.init_logger(self.base_path)
        self.load(config_path)

    def register_folder(self, name):
        # Separate the naming and construction of folders
        new_path = os.path.join(self.output_path, name)
        self._output_folders.append(new_path)
        setattr(self, name + "_path", new_path)

    def make_output_folders(self):
        util.make_dir(self.output_path)
        for pth in self._output_folders:
            util.make_dir(pth)

    def register_output_folders(self):
        self.register_folder('subsets')
        self.register_folder('schemes')
        self.register_folder('phylofiles')
        self.register_folder('start_tree')

    def init_logger(self, pth):
        handler = logging.FileHandler(os.path.join(pth, "log.txt"), 'a')
        formatter = logging.Formatter(
            "%(levelname)-8s | %(asctime)s | %(name)-10s | %(message)s")
        handler.setFormatter(formatter)
        handler.setLevel(logging.DEBUG)
        logging.getLogger("").addHandler(handler)

    def load(self, config_path):
        """We get the parser to construct the configuration"""
        log.info("Loading configuration at '%s'", config_path)
        self.config_path = config_path
        p = parser.Parser(self)
        p.parse_file(config_path)
        log.info(
            "------------------------ BEGINNING NEW RUN -------------------------------")

    def set_base_path(self, base_path):
        self.full_base_path = os.path.abspath(base_path)
        log.info("Setting working folder to: '%s'", self.full_base_path)

        # Now make our working folder this folder. All of our other paths will
        # be relative to this
        self.old_working_directory = os.getcwd()
        os.chdir(self.full_base_path)

        # Our base path is now this
        self.base_path = '.'
        self.output_path = os.path.join(self.base_path, "analysis")
        self.full_output_path = os.path.join(self.full_base_path, "analysis")

    def set_alignment_file(self, align):
        log.info("Setting 'alignment' to '%s'", align)
        self.alignment = align

    def set_option(self, option, value):
        # Make everything lowercase, this makes life easier for us
        value = value.lower()

        if option not in self.options:
            log.error("'%s' is not a valid option to set in the configuration",
                      option)
            raise ConfigurationError

        # Compare lower case
        valid = [x.lower() for x in self.options[option]]
        if value not in valid:
            log.error("'%s' is not a valid option for '%s'" % (value, option))
            log.info("The only valid options for '%s' are: %s" %
                     (option, "'%s'" % ("', '".join(self.options[option]))))
            raise ConfigurationError

        # TODO: not the best place for this at all..., but it works
        if option == "search" and "cluster" in value and self.phylogeny_program != 'raxml':
            log.error("Clustering methods are only available when using raxml"
                      " (the --raxml commandline option). Please check and try again."
                      " See the manual for more details.")
            raise ConfigurationError

        if option == "search" and "kmeans" in value and self.phylogeny_program != 'phyml' and self.kmeans_opt != 1:
            log.error("You have chosen a kmeans option (--kmeans_opt) that does not work "
                      "with the --raxml option. Please re-run your analysis in one of two ways: "
                      "\n 1. Remove the --raxml commandline option, so that PhyML is used, or "
                      "\n 2. Change the --kmeans_opt commandline option to 1 (or remove it) and leave the "
                      "--raxml option in place.")
            raise ConfigurationError

        log.info("Setting '%s' to '%s'", option, value)
        setattr(self, option, value)

    def validate(self):
        """Should be called before processing"""
        # Just path validation for now.
        util.check_folder_exists(self.base_path)
        self.alignment_path = os.path.join(self.base_path, self.alignment)
        log.info("Looking for alignment file '%s'...", self.alignment_path)
        util.check_file_exists(self.alignment_path)
        if self.user_tree is None:
            self.user_tree_topology_path = None
        else:
            self.user_tree_topology_path = \
                    os.path.join(self.base_path, self.user_tree)
            log.info("Looking for tree file '{}'...".format(
                self.user_tree_topology_path))
            util.check_file_exists(self.user_tree_topology_path)

