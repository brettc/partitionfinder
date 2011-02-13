import logging
log = logging.getLogger("config")

import os

from partition import Partition
from parser import Parser, ParserError

from pyparsing import (
    Word, Dict, OneOrMore, alphas, nums, Suppress, Optional, Group, stringEnd,
    delimitedList, pythonStyleComment, ParseException, line, lineno, col,
    Keyword)

class ConfigurationError(Exception):
    pass

class Configuration(object):
    """We use this to hold all of the configuration info"""
    def __init__(self, base_path):
        self.base_path = os.path.abspath(base_path)
        if not os.path.exists(self.base_path) or \
           not os.path.isdir(self.base_path):
            log.error("No such folder: '%s'", self.base_path)
            raise ConfigurationError

        self.config_path = os.path.join(self.base_path, "partition_finder.cfg")
        self.output_path = os.path.join(self.base_path, "output")
        self.log_path = os.path.join(self.base_path, "partition_finder.log")

        # Add a log file in this folder
        log.info("Using folder: '%s'", self.base_path)
        self.init_log()

    def load(self):
        if not os.path.exists(self.config_path) or \
           not os.path.isfile(self.config_path):
            log.error("Configuration file '%s' does not exist",
                      self.config_path)
            raise ConfigurationError

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







