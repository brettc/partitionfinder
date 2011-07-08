import logging
log = logging.getLogger("config")

import os

import scheme, subset, partition
import util
class ConfigurationError(util.PartitionFinderError):
    pass

class Configuration(object):
    """This holds the user configuration info"""

    # List of valid options. The first one is the default
    options = {
        'branchlengths': ['linked', 'unlinked'],
        'model_selection': ['AIC', 'AICc', 'BIC'],
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
        if value not in valid:
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

if __name__ == '__main__':
    import logging
    logging.basicConfig(level=logging.DEBUG)
    # report_settings()
    
