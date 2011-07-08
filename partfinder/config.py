import logging
log = logging.getLogger("config")

import os

import scheme, subset, partition
import util
class ConfigurationError(util.PartitionFinderError):
    pass

class Configuration(object):
    """This holds the user configuration info"""
    def __init__(self):
        self.partitions = partition.PartitionSet()
        self.schemes = scheme.SchemeSet()
        
    def set_base_path(self, base_path):
        log.info("Using folder: '%s'", base_path)
        self.base_path = base_path
        self.output_path = os.path.join(base_path, "analysis")

    def set_alignment_file(self, align):
        log.info("Setting 'alignment' to '%s'", align)
        self.alignment = align
        self.alignment_path = os.path.join(self.base_path, align)

    def validate(self):
        """Should be called before processing"""
        util.check_folder_exists(self.base_path)
        util.check_file_exists(self.alignment_path)
        # settings.alignment_path = os.path.join(settings.base_path,
                                            # settings.alignment)
        # _check_file(settings.alignment_path)

    # def __getattr__(self, name):
        # if name not in self.__dict__:
            # log.error("The setting '%s' is not defined. "
                      # "Did you initialise the configuration?", 
                      # name)
            # raise ConfigurationError

        # return self.__dict__[name]

# def report_settings():
    # log.debug("Settings are as follows:")
    # for x in settings.__dict__:
        # if not x.startswith('__'):
            # log.debug("%s: %s", x, getattr(settings, x))

if __name__ == '__main__':
    import logging
    logging.basicConfig(level=logging.DEBUG)
    # report_settings()
    
