import logging
log = logging.getLogger("config")

import os, shutil

from partition import Partition
from parser import Parser, ParserError
from fasta import read_fasta, FastaError

__all__ =  ["ConfigurationError", "settings", "initialise"]

class ConfigurationError(Exception):
    pass

class ProcessingError(Exception):
    pass

def _check_file(pth):
    if not os.path.exists(pth) or not os.path.isfile(pth):
        log.error("No such file: '%s'", pth)
        raise ConfigurationError

def _check_folder(pth):
    if not os.path.exists(pth) or not os.path.isdir(pth):
        log.error("No such folder: '%s'", pth)
        raise ConfigurationError

def _make_folder(pth):
    if os.path.exists(pth):
        if not os.path.isdir(pth):
            log.error("Cannot create folder '%s'", pth)
            raise ConfigurationError
    else:
        os.mkdir(pth)


class Settings(object):
    """This holds the configuration info"""
    def __init__(self):
        pass

    def __getattr__(self, name):
        if name not in self.__dict__:
            log.error("The setting '%s' is not defined. "
                      "Did you initialise the configuration?", 
                      name)
            raise ConfigurationError

        return self.__dict__[name]

init_done = False
settings = Settings()

def initialise(pth, force_restart=False):
    global init_done
    if init_done:
        log.error("Cannot initialise more than once")
        raise ConfigurationError

    # Allow for user and environment variables
    pth = os.path.expanduser(pth)
    pth = os.path.expandvars(pth)
    pth = os.path.normpath(pth)
    settings.base_path = pth
    settings.force_restart = force_restart

    _check_folder(settings.base_path)
    log.info("Using folder: '%s'", settings.base_path)

    create_debug_log()

    settings.output_path = os.path.join(settings.base_path, "output")
    if settings.force_restart:
        if os.path.exists(settings.output_path):
            log.warning("Deleting all previous workings in '%s'", settings.output_path)
            shutil.rmtree(settings.output_path)

    _make_folder(settings.output_path)

    find_modelgenerator()

    init_done = True

def create_debug_log():
    """Add a full debug log file in the folder"""
    settings.log_path = os.path.join(settings.base_path, "partition_finder.log")
    log.info("Full log is in: '%s'", settings.log_path)

    # Append to the log file. we'll get multiple runs then
    log_output = logging.FileHandler(settings.log_path, mode='a')
    log_output.setLevel(logging.DEBUG)
    formatter = logging.Formatter(
        '%(levelname)-8s | %(asctime)s | %(name)-10s | %(message)s',
        datefmt="%Y-%m-%d %H:%M:%S")
    log_output.setFormatter(formatter)
    logging.getLogger('').addHandler(log_output)
    # Mark a new session so it is easy to read in the log file
    log.debug("------------------ NEW LOG SESSION BEGINS -----------------")

def load():
    settings.config_path = os.path.join(settings.base_path, "partition_finder.cfg")
    _check_file(settings.config_path)

    log.info("Loading configuration at '%s'", settings.config_path)
    try:
        p = Parser(self)
        self.processing = p.parse_file(self.config_path)
    except ParserError, p:
        # Catch any parsing errors, print something out, then raise a
        # configuration error, as this is the general error we expect from
        # this part of the process
        log.error(p.format_message())
        raise ConfigurationError

    self.alignment_path = os.path.join(self.base_path,
                                        self.alignment_file)
    _check_file(self.alignment_path)

    # Now read in the sequence as part of the loading
    try:
        self.sequence = read_fasta(self.alignment_path)
    except FastaError:
        log.error("Cannot load Fasta file '%s'" % self.alignment_path)
        raise ConfigurationError

def find_modelgenerator():
    """Make sure we know where the java file is..."""

    pth = os.path.abspath(__file__)
    # Split off the name and the directory...
    pth, notused = os.path.split(pth)
    pth, notused = os.path.split(pth)
    # Now go back down into programs...
    pth = os.path.join(pth, "programs", "modelgenerator.jar")
    pth = os.path.normpath(pth)

    log.debug("Checking for modelgenerator program")
    _check_file(pth)
    log.debug("Modelgenerator program found at '%s'" % pth)
    settings.modelgen_path = pth

if __name__ == '__main__':
    import logging
    logging.basicConfig(level=logging.DEBUG)
    initialise("~/tmp", True)
