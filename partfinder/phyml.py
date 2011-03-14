"""Run phyml and parse the output"""

import logging
log = logging.getLogger("phyml")

import config

import subprocess, shlex, tempfile, os, shutil, sys

from pyparsing import (
    Word, Literal, alphas, nums, Suppress, Group, stringEnd, ParseException,
    line, lineno, col, LineStart, SkipTo, LineEnd,
    )

from phyml_models import get_model_commandline

class PhymlError(Exception):
    pass

def run_phyml(command):
    log.debug("Running command '%s'", command)

    # Note: We use shlex.split as it does a proper job of handling command
    # lines that are complex
    try:
        subprocess.check_call(shlex.split(command), shell=False)
	
    except subprocess.CalledProcessError:
        log.error("command '%s' failed to execute successfully", command)
        raise PhymlError

def dupfile(src, dst):
    # Make a copy or a symlink so that we don't overwrite different model runs
    # of the same alignment
    
    # TODO maybe this should throw...?
    if os.path.exists(dst):
        os.remove(dst)
    if sys.platform in ("darwin", "unix"):
        os.symlink(src, dst)
    else:
        shutil.copyfile(src, dst)

def mvfile(src, dst):
    pass

# This should generate a bootstrap tree 
def make_tree(program, alignment):
    # We can get a rough (but probably correct!) topology using BioNJ, this is like
    # Neighbour joining but a little bit better. Turns out we can do this already
    # with PhyML like this:
    log.debug("Making a tree for %s", alignment)

    # TODO: put into a temp folder?

    # First get the BioNJ topology like this:
    command = "%s -i %s -o n -b 0" % (program, alignment)
    run_phyml(command)
    output_path = tree_path(alignment)

    # in our case (well behaved data) we get the right topology already like
    # this, and nearly correct branchlengths too. Even so, we might want to
    # re-estimate branchlengths from scratch using a better model.

    # To do that, first rename the tree 
    dirpth, fname = os.path.split(output_path)
    tree_pth = os.path.join(dirpth, 'bionj_tree.phy')
    log.debug("Moving %s to %s", output_path, tree_pth)
    os.rename(output_path, tree_pth) 

    # now, we run PhyML asking it to re-optimise the branch lengths according to
    # some model, but not to mess with the topology.
    #
    # using "-o lr" will optimise the model parameters and brlens together. This
    # is important since the two interact, i.e. good brlens require good model
    # parameters.
    command = "%s -i %s -u %s -m GTR -c 4 -a e -v e -o lr -b 0" % (
        program, alignment, tree_pth)
    run_phyml(command)

    # Now return the path of the final tree alignment
    return tree_path(alignment)


def analyse(program, model, alignment, tree):
    """Do the analysis -- this will overwrite stuff!"""

    analysis, output = analysis_path(alignment, model)

    dupfile(alignment, analysis)

    model_params = get_model_commandline(model)

    command = "%s -i %s -u %s %s" % (program, analysis, tree, model_params)
    run_phyml(command)

    # Now get rid of this -- we have the original 
    os.remove(analysis)

    # Return the output path, where we can read the info
    return output

def tree_path(alignment):
    pth, ext = os.path.splitext(alignment)
    return pth + ".phy_phyml_tree.txt"

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
        return "PhymlResult(model:%s, lnl:%s, secs:%s)" % (
            self.model, self.lnl, self.seconds) 

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

# Stateless, so safe
the_parser = Parser()

def parse(text):
    return the_parser.parse(text)

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)

    import tempfile

    tmp = tempfile.mkdtemp()
    config.initialise(tmp)

    s = config.settings
    al = os.path.join(s.test_path, 'part1.phy')
    make_tree(s.program_path, al)


def test1():
    import phyml_models 
    import tempfile

    tmp = tempfile.mkdtemp()
    config.initialise(tmp)

    s = config.settings
    # TODO Copy to temp
    al = os.path.join(s.test_path, 'part1.phy')

    tr = os.path.join(s.test_path, 'start_tree_correct.phy')

    for model in phyml_models.get_all_models():
        print "Analysing using model %s:" % model,
        out_pth = analyse(s.program_path, model, al, tr)
        output = open(out_pth, 'rb').read()
        res = parse(output)
        print res.lnl


    shutil.rmdir(tmp)
