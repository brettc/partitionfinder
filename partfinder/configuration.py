import logging
log = logging.getLogger("config")

import os

from partition import Partition
from parser import Parser, ParserError
from fasta import read_fasta, FastaError

from pyparsing import (
    Word, Dict, OneOrMore, alphas, nums, Suppress, Optional, Group, stringEnd,
    delimitedList, pythonStyleComment, ParseException, line, lineno, col,
    Keyword)

class ConfigurationError(Exception):
    pass

def _check_file(pth):
    if not os.path.exists(pth) or not os.path.isfile(pth):
        log.error("No such file: '%s'", pth)
        raise ConfigurationError

def _check_folder(pth):
    if not os.path.exists(pth) or not os.path.isdir(pth):
        log.error("No such folder: '%s'", pth)
        raise ConfigurationError

class Configuration(object):
    """We use this to hold all of the configuration info"""
    def __init__(self, base_path):
        self.base_path = os.path.abspath(base_path)
        _check_folder(self.base_path)

        self.config_path = os.path.join(self.base_path, "partition_finder.cfg")
        self.output_path = os.path.join(self.base_path, "output")
        self.log_path = os.path.join(self.base_path, "partition_finder.log")

        # Add a log file in this folder
        log.info("Using folder: '%s'", self.base_path)
        self.init_log()

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
        log.info("Creating full log in: '%s'", self.log_path)
        log_output = logging.FileHandler(self.log_path, mode='w')
        log_output.setLevel(logging.DEBUG)
        formatter = logging.Formatter(
            '%(levelname)-8s | %(asctime)s | %(name)-10s | %(message)s',
            datefmt="%Y-%m-%d %H:%M:%S")
        log_output.setFormatter(formatter)
        logging.getLogger('').addHandler(log_output)
        # logging.getLogger('').setLevel(logging.DEBUG)
        log.debug("Log file created at '%s'", self.log_path)

    def has_changed(self):
        """check to see if the config has changed from when it was
        last run"""
        pass
        # Write a signature out somewhere...







