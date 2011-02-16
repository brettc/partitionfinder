import logging
log = logging.getLogger("config")

import os, shutil

from partition import Partition
from parser import Parser, ParserError
from fasta import read_fasta, FastaError

from pyparsing import (
    Word, Dict, OneOrMore, alphas, nums, Suppress, Optional, Group, stringEnd,
    delimitedList, pythonStyleComment, ParseException, line, lineno, col,
    Keyword)

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

class Configuration(object):
    """We use this to hold all of the configuration info"""
    def __init__(self, base_path, force_restart):
        self.base_path = os.path.abspath(base_path)
        self.force_restart = force_restart
        _check_folder(self.base_path)
        # Add a log file in this folder
        log.info("Using folder: '%s'", self.base_path)

        self.log_path = os.path.join(self.base_path, "partition_finder.log")
        self.init_log()

        # Check all the folders and files we need...
        self.config_path = os.path.join(self.base_path, "partition_finder.cfg")

        self.output_path = os.path.join(self.base_path, "output")
        if self.force_restart:
            log.warning("Deleting all previous workings in '%s'", self.output_path)
            shutil.rmtree(self.output_path)

        _make_folder(self.output_path)

        self.find_modelgenerator()

    def find_modelgenerator(self):
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

        self.modelgen_path = pth


    def load(self):
        _check_file(self.config_path)

        log.info("Loading configuration at '%s'", self.config_path)
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

    def init_log(self):
        """Add a full debug log file in the folder"""
        log.info("Using full log in: '%s'", self.log_path)

        # Append to the log file. we'll get multiple runs then
        log_output = logging.FileHandler(self.log_path, mode='a')
        log_output.setLevel(logging.DEBUG)
        formatter = logging.Formatter(
            '%(levelname)-8s | %(asctime)s | %(name)-10s | %(message)s',
            datefmt="%Y-%m-%d %H:%M:%S")
        log_output.setFormatter(formatter)
        logging.getLogger('').addHandler(log_output)

    def has_changed(self):
        """Check if changed from previous run
        """
        # TODO Write a signature out somewhere...

    def process(self):
        # TODO: This is one way to do it...
        # Or we could do it by scheme?
        
        # Go through all the subsets and run them
        for ss_id, ss in self.schemes.subsets.iteritems():
            ss.process(self)

        # for sch in self.schemes:
            # scheme.assess()



