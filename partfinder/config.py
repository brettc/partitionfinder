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
log = logging.getLogger("config")

import os
import fnmatch
import cPickle as pickle
import scheme
import partition
import parser
import util
import progress
import numbers


class ConfigurationError(util.PartitionFinderError):
    pass


class Configuration(object):
    """This holds the user configuration info"""

    # List of valid options. The first one is the default
    options = {
        'branchlengths': ['linked', 'unlinked'],
        'model_selection': ['aic', 'aicc', 'bic'],
        'search': ['all', 'user', 'greedy', 'clustering', 'greediest']
    }

    def __init__(self, datatype="DNA", phylogeny_program='phyml',
        save_phylofiles=False, cmdline_extras = "", cluster_weights = None, 
        greediest_schemes = 1, greediest_percent=1):

        self.partitions = partition.PartitionSet()
        self.schemes = scheme.SchemeSet()
        self.save_phylofiles = save_phylofiles
        self.progress = progress.NoProgress(self)
        self.cmdline_extras = cmdline_extras
        self.greediest_schemes = greediest_schemes
        self.greediest_percent = float(greediest_percent)

        self.base_path = '.'
        self.alignment = None
        self.user_tree = None

        ### Some basic checking of the setup, so that we don't hit too many problems later
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

        if cluster_weights==None:
            #default to equal weights. TODO. This should change depending on results of our analyses
            self.cluster_weights = {"rate": 1, "freqs": 1, "model": 1, "alpha": 1}
        else:
            #TODO. Is there a more robust way to do this...
            cluster_weights = cluster_weights.split(",")

            #now we check that it's a list of exactly three numbers
            if len(cluster_weights)!=4:
                log.error("Your --cluster_weights argument should have exactly 4"
                    " numbers separated by commas, but it has %d. "
                    "Please check and try again" %len(cluster_weights))
                raise ConfigurationError
            final_weights = []
            for thing in cluster_weights:
                try:
                    num = float(eval(thing))
                    final_weights.append(num)
                except:
                    log.error("Unable to understand your --cluster_weights argument."
                        " It should look like this: --cluster_weights '1,2,3,6'. "
                        "Please double check that you included quotes, "
                        "and three numbers separated by commas. Then try again. "
                        "The part that I coudln't understand is this: '%s'" %thing)
                    raise ConfigurationError

            log.info("Setting cluster_weights to: subset_rate = %.1f, freqs = %.1f, model = %.1f" %(final_weights[0], final_weights[1], final_weights[2]))
            self.cluster_weights = {}
            self.cluster_weights["rate"] = final_weights[0]
            self.cluster_weights["freqs"] = final_weights[1]
            self.cluster_weights["model"] = final_weights[2]
            self.cluster_weights["alpha"] = final_weights[3]


        # Set the defaults into the class. These can be reset by calling
        # set_option(...)
        for o, v in self.options.items():
            # Could call self.set_option here -- but it might confuse users
            setattr(self, o, v[0])

    def reset(self):
        # Only required in tests, when we keep changing folders
        if hasattr(self, 'old_cwd'):
            os.chdir(self.old_cwd)

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

        log.error("Cannot find a configuration file"
                  "working folder '%s'", folder)

        raise ConfigurationError

    def load_base_path(self, pth):
        """Load using a base path folder"""
        # Allow for user and environment variables
        pth = os.path.expanduser(pth)
        pth = os.path.expandvars(pth)
        pth = os.path.normpath(pth)

        folder, filename = self.find_config_file(pth)

        #check that user didn't enter a file instead of a folder
        # if os.path.isfile(pth):
            # log.error("The second argument of the commandline currently points to a file, but it should point to the folder that contains the alignment and .cfg files, please check.")
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

        self.make_output_folders()

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
        self.register_folder('phyml')
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
        log.info("------------------------ BEGINNING NEW RUN -------------------------------")

    def set_base_path(self, base_path):
        self.full_base_path = os.path.abspath(base_path)
        log.info("Setting working folder to: '%s'", self.full_base_path)

        # Now make our working folder this folder. All of our other paths will
        # be relative to this
        self.old_cwd = os.getcwd()
        os.chdir(base_path)

        # Our base path is now this
        self.base_path = '.'
        self.output_path = os.path.join(self.base_path, "analysis")
        self.full_output_path = os.path.join(self.full_base_path, "analysis")

    def set_alignment_file(self, align):
        log.info("Setting 'alignment' to '%s'", align)
        self.alignment = align

    def set_option(self, option, value):

        #make everything lowercase, this makes life easier for us
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

        #TODO: not the best place for this at all..., but it works
        if option=="search" and value=="clustering" and self.phylogeny_program!='raxml':
            log.error("The 'search = clustering' option is only availalbe when using raxml"
                " (the --raxml commandline option). Please check and try again."
                " See the manual for more details.")
            raise ConfigurationError

        #TODO: not the best place for this at all..., but it works
        if option=="search" and value=="greediest" and self.phylogeny_program!='raxml':
            log.error("The 'search = greediest' option is only availalbe when using raxml"
                " (the --raxml commandline option). Please check and try again."
                " See the manual for more details.")
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
            self.user_tree_topology_path = os.path.join(self.base_path,
                                                        self.user_tree)
            log.info(
                "Looking for tree file '%s'...", self.user_tree_topology_path)
            util.check_file_exists(self.user_tree_topology_path)

    def check_for_old_config(self):
        """Check whether the analysis dictated by cfg has been run before, and if the config has changed
        in any way that would make re-running it invalid"""
        #the important stuff in our analysis, that can't change if we want to re-use old subsets
        if self.user_tree is None:
            topology = ""
        else:
            topology = open(self.user_tree_topology_path).read()

        cfg_list = [self.alignment,
                    self.branchlengths,
                    self.partitions.partitions,
                    self.phylogeny_program,
                    topology]

        #we need to know if there's anything in the subsets folder
        subset_path = os.path.join(self.output_path, 'subsets')
        has_subsets = False
        if os.path.exists(subset_path):
            for file in os.listdir(subset_path):
                if fnmatch.fnmatch(file, '*.bin'):
                    has_subsets = True
                    break

        #we also need to know if there's an old conifg file saved
        cfg_dir = os.path.join(self.output_path, 'cfg')
        old_cfg_path = os.path.join(cfg_dir, 'oldcfg.bin')
        if os.path.exists(old_cfg_path):
            has_config = True
        else:
            has_config = False

        if has_subsets == False:
            #we have no subsets so can't screw anything up, just copy the new cfg file settings, overwrite anything else
            if not os.path.exists(cfg_dir):
                os.makedirs(cfg_dir)
            #store a nice binary
            f = open(old_cfg_path, 'wb')
            pickle.dump(cfg_list, f, -1)
            f.close()
            return 0

        else:  # there are subsets
            if has_config == False:
                log.error("There are subsets stored, but PartitionFinder can't determine where they are from")
                log.info("Please re-run the analysis using the '--force-restart' option at the command line")
                log.warning("This will delete all of the analyses in the '/analysis' folder")
                raise ConfigurationError
            else:
                #we have an old config, load it and compare the important bits
                f = open(old_cfg_path, 'rb')
                old_cfg = pickle.load(f)
                f.close()
                fail = []

                if len(old_cfg) != len(cfg_list):
                    log.error("Your old configuration doesn't match with your new one")
                    log.error("The most common cause of this error is trying to run a half-finished analysis"
                              " but switching PartitionFinder versions half way through")
                    log.error("The solution is to either go back to the same version of PartitionFinder"
                              " you used for the initial analysis, or to re-run the analysis using the "
                              "--force-restart option at the command line. Note that this will delete "
                              "all previous analyses in the '/analysis' folder")

                if not old_cfg[0] == cfg_list[0]:
                    fail.append("alignment")
                if not old_cfg[1] == cfg_list[1]:
                    fail.append("branchlengths")

                old_parts = set([str(part) for part in old_cfg[2]])
                new_parts = set([str(part) for part in cfg_list[2]])
                if len(old_parts.difference(new_parts)) > 0:
                    fail.append("[data_blocks]")

                if not old_cfg[3] == cfg_list[3]:
                    fail.append(
                        "phylogeny_program (the --raxml commandline option)")

                if not old_cfg[4] == cfg_list[4]:
                    fail.append("user_tree_topology")

                if len(fail) > 0:
                    log.error("There are subsets stored, but PartitionFinder has detected that these were run using a different .cfg setup")
                    log.error("The following settings in the new .cfg file are incompatible with the previous analysis: %s" % (', '.join(fail)))
                    log.info("To run this analysis and overwrite previous output, re-run the analysis using '--force-restart' option")
                    log.info("To run this analysis without deleting the previous analysis, please place your alignment and .cfg in a new folder and try again")
                    raise ConfigurationError
