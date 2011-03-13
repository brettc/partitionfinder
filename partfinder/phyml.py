"""Run phyml and parse the output"""

import logging
log = logging.getLogger("phyml")

import config

import subprocess, shlex, tempfile, os, shutil, sys

from pyparsing import (
    Word, Literal, alphas, nums, Suppress, Group, stringEnd, ParseException,
    line, lineno, col, LineStart, SkipTo, LineEnd,
    )

class PhymlError(Exception):
    pass

# This should generate a bootstrap tree 
def make_tree():
    pass

def analyse(program, model, alignment, tree):

    analysis, output = analysis_path(alignment, model)

    # Make a copy or a symlink so that we don't overwrite different model runs
    # of the same alignment
    if sys.platform in ("darwin", "unix"):
        os.symlink(alignment, analysis)
    else:
        shutil.copyfile(alignment, analysis)

    command = "%s -i %s -u %s -m %s -c 8 -a e -o lr -b 0 --constrained_lens" % (
        program, analysis, tree, model)

    log.debug("Running command '%s'", command)

    # Note: We use shlex.split as it does a proper job of handling command
    # lines that are complex
    p = subprocess.Popen(
        shlex.split(command),
        shell=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)

    # Capture the output, we might put it into the errors
    stdout, stderr = p.communicate()

    if p.returncode != 0:
        log.error("program failed to execute successfully: output follows")
        log.error(stderr)
        raise PhymlError

    # Now get rid of this -- we have the original 
    os.remove(analysis)

    # Return the output path, where we can read the info
    return output

def analysis_path(alignment, model):
    pth, ext = os.path.splitext(alignment)

    analyse = pth + "[" + model + "]" + ext

    pth, ext = os.path.splitext(analyse)
    output = pth + ".phy_phyml_stats.txt"

    return analyse, output

class PhymlResult(object):
    def __init__(self, model, lnl, seconds):
        self.lnl = lnl
        self.model = model
        self.seconds = seconds

    def __str__(self):
        return "PhymlResult(model:%s, lnl:%s, secs:%s)" % (self.model,
                                                           self.lnl,
                                                           self.seconds)

class Parser(object):
    def __init__(self):
        FLOAT = Word(nums + '.-').setParseAction(lambda x: float(x[0]))
        INTEGER = Word(nums + '-').setParseAction(lambda x: int(x[0]))

        OB = Suppress("(")
        CB = Suppress(")")
        LNL_LABEL = Literal("Log-likelihood:")
        TIME_LABEL = Literal("Time used:")
        HMS = Word(nums + "hms") # A bit rough...
        MODEL_LABEL = Literal("Model of nucleotides substitution:")
        MODEL = Word(alphas + nums)

        lnl = (LNL_LABEL + FLOAT("lnl"))
        time = (TIME_LABEL + HMS("time") + OB + INTEGER("seconds") + Suppress("seconds") + CB)
        model = (MODEL_LABEL + MODEL("model"))

        # Shorthand...
        def nextbit(label, val):
            return Suppress(SkipTo(label)) + val

        # Just look for these things
        self.root_parser = \
                nextbit(MODEL_LABEL, model) +\
                nextbit(LNL_LABEL, lnl) +\
                nextbit(TIME_LABEL, time)

    def parse(self, text):
        log.info("Parsing phyml output...")
        try:
            tokens = self.root_parser.parseString(text)
        except ParseException, p:
            log.error(p.format_message())
            raise PhymlError

        result = PhymlResult(
            lnl=tokens.lnl, seconds=tokens.seconds, model=tokens.model)

        log.info("Parsed phyml output is %s", result)
        return result

        # We'll scan for what we want, rather than parse the whole thing
        # for match in self.match_parser.scanString(text):
            # tokens = match[0]
            # # print tokens
            # if tokens.lnl: print tokens.lnl
            # if tokens.time: print tokens.time[2]

# Stateless, so safe
the_parser = Parser()

def parse(text):
    return the_parser.parse(text)

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    config.initialise("~/tmp")
    s = config.settings
    al = os.path.join(s.test_path, 'part1.phy')
    tr = os.path.join(s.test_path, 'start_tree_correct.phy')
    out_pth = analyse(s.program_path, "HKY", al, tr)
    output = open(out_pth, 'rb').read()
    p = Parser()
    res = p.parse(output)
