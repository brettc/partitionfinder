import logging
log = logging.getLogger("modelgen")

import config
import subprocess

from pyparsing import (
    Word, Dict, OneOrMore, alphas, nums, Suppress, Optional, Group, stringEnd,
    delimitedList, pythonStyleComment, ParseException, line, lineno, col,
    Keyword, ParserElement, ParseException)

def run_modelgen(alignment, num_gamma_cat=4):
    mg = config.settings.modelgen_path
    command = "java -jar %s %s %d" % (mg, alignment_path, num_gamma_cat)
    p = subprocess.Popen(command.split(),
                         shell=False,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    import config
    config.initialise("~/tmp")
    run_modelgen("xxxx")

