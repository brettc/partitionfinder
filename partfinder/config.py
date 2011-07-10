import logging
log = logging.getLogger("config")

import os, fnmatch
import cPickle as pickle
import scheme, subset, partition
import util
class ConfigurationError(util.PartitionFinderError):
    pass

class Configuration(object):
    """This holds the user configuration info"""

    # List of valid options. The first one is the default
    options = {
        'branchlengths': ['linked', 'unlinked'],
        'model_selection': ['aic', 'aicc', 'bic'],
        'search': ['all', 'user', 'greedy']
        }

    def __init__(self):
        self.partitions = partition.PartitionSet()
        self.schemes = scheme.SchemeSet()

        self.base_path = '.'
        self.alignment = None

        # Set the defaults into the class. These can be reset by calling
        # set_option(...)
        for o, v in self.options.items():
            # Could call self.set_option here -- but it might confuse users
            setattr(self, o, v[0])
            
    def set_base_path(self, base_path):
        log.info("Using folder: '%s'", base_path)
        self.base_path = base_path
        self.output_path = os.path.join(base_path, "analysis")

    def set_alignment_file(self, align):
        log.info("Setting 'alignment' to '%s'", align)
        self.alignment = align

    def set_option(self, option, value):
        if option not in self.options:
            log.error("'%s' is not a valid option to set in the configuration",
                      option)
            raise ConfigurationError

        valid = self.options[option]
     
        #force lowercase for comparisons
        value = value.lower()
        valid_lower = []
        [valid_lower.append(thing.lower()) for thing in valid]        
        if value not in valid_lower:
            log.error("'%s' is not a valid option for '%s'" % (value, option))
            log.info("The only valid options for '%s' are: %s" % 
                     (option, "'%s'" %("', '".join(self.options[option]))))
            raise ConfigurationError

        log.info("Setting '%s' to '%s'", option, value)
        setattr(self, option, value)

    def validate(self):
        """Should be called before processing"""
        # Just path validation for now.
        # TODO: test the alignment against the ranges in the partitions
        util.check_folder_exists(self.base_path)
        self.alignment_path = os.path.join(self.base_path, self.alignment)
        util.check_file_exists(self.alignment_path)

def check_for_old_config(cfg):
    """Check whether the analysis dictated by cfg has been run before, and if the config has changed
    in any way that would make re-running it invalid"""
    #the important stuff in our analysis, that can't change if we want to re-use old subsets
    cfg_list = [cfg.alignment, cfg.branchlengths, cfg.partitions.partitions, cfg.models]

    #we need to know if there's anything in the subsets folder
    subset_path = "%s/subsets/" %(cfg.output_path)
    has_subsets = False
    if os.path.exists(subset_path):
        for file in os.listdir(subset_path):
            if fnmatch.fnmatch(file, '*.bin'):
                has_subsets=True
                break

    #we also need to know if there's an old conifg file saved
    cfg_dir = "%s" %(cfg.output_path) 
    old_cfg_path = "%s/oldcfg.bin" %(cfg_dir)
    if os.path.exists(old_cfg_path):
        has_config = True
    else:
        has_config = False

    if has_subsets==False:
        #we have no subsets so can't screw anything up, just copy the new cfg file settings, overwrite anything else
        if not os.path.exists(cfg_dir):
            os.mkdir(cfg_dir)
        #store a nice binary
        f = open(old_cfg_path, 'wb')
        pickle.dump(cfg_list, f, -1)
        return 0

    else: #there are subsets
        if has_config==False:
            log.error("There are subsets stored, but PartitionFinder can't determine where they are from")
            log.info("Please re-run the analysis using the '--force-restart' option at the command line")
            log.warning("This will delete all of the analyses in the '/analysis' folder")
            raise ConfigurationError
        else:
            #we have an old config, load it and compare the important bits
            f = open(old_cfg_path, 'rb')
            old_cfg = pickle.load(f)
            fail = []
            
            if not old_cfg[0]==cfg_list[0]:
                fail.append("alignment")
            if not old_cfg[1]==cfg_list[1]:
                fail.append("branchlengths")
            if not old_cfg[3]==cfg_list[3]:
                fail.append("models")
            
            old_parts = set()
            [old_parts.add(str(part)) for part in old_cfg[2]]
            new_parts = set()
            [new_parts.add(str(part)) for part in cfg_list[2]]
            if len(old_parts.difference(new_parts))>0:
                fail.append("[partitions]")
            
            if len(fail)>0:
                log.error("There are subsets stored, but PartitionFinder has detected that these were run using a different .cfg setup")
                log.error("The following settings in the new .cfg file are incompatible with the previous analysis: %s" %(', '.join(fail)))
                log.info("To run using the current .cfg file, and overwrite previous output, please re-run the analysis using the '--force-restart' option at the command line")
                log.info("To run a different analysis, please place your alignment and .cfg in a new folder and try again")
                raise ConfigurationError



if __name__ == '__main__':
    import logging
    logging.basicConfig(level=logging.DEBUG)
    # report_settings()
    
