"""Run phyml and parse the output"""

import logging
log = logging.getLogger("phyml")

import config

import subprocess, shlex, tempfile, os, shutil

from pyparsing import (
    Word, Literal, alphas, nums, Suppress, Group, stringEnd, ParseException,
    line, lineno, col, LineStart, SkipTo, LineEnd,
    )

class PhymlError(Exception):
    pass

def run(params):
    # Run it in a temporary directory, so we can ditch all the output.
    output_dir = tempfile.mkdtemp()
    log.debug("Created temporary folder %s for phyml output", output_dir)

    model = "HKY"
    binary_path = config.settings.program_path
    alignment_path = "xxx.phy"
    tree_path = "xxx_tree.phy"
    # How does this work?
    output_path = alignment + ".phy_phyml_stats.txt"

    command = "%s -i %s -u %s -m %s -c 8 -a e  --free_rates -o r -b 0" % (
        binary_path, alignment_path, tree_path, model)

    try:
        log.debug("Running command '%s'", command)

        # Note: We use shlex.split as it does a proper job of handling command
        # lines that are complex
        p = subprocess.Popen(
            shlex.split(command),
            shell=False,
            cwd=output_dir,
            stdout=subprocess.PIPE)
        stdout, stderr = p.communicate()
        output = open(output_path, 'rb').read()
    finally:
        log.debug("Removing temporary folder %s", output_dir)
        shutil.rmtree(output_dir)

    if p.returncode != 0:
        log.error("Modelgenerator program failed to execute successfully")
        raise ModelGeneratorError

    # Otherwise this should contain the output. It will need parsing
    return output

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
        def NextLabel(label, val):
            return Suppress(SkipTo(label)) + val

        # Just look for these things
        self.root_parser = \
                NextLabel(MODEL_LABEL, model) +\
                NextLabel(LNL_LABEL, lnl) +\
                NextLabel(TIME_LABEL, time)

    def parse(self, text):
        try:
            tokens = self.root_parser.parseString(text)
        except ParseException, p:
            log.error(p.format_message())
            raise PhymlError

        return PhymlResult(
            lnl=tokens.lnl, 
            seconds=tokens.seconds, 
            model=tokens.model)

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

typical_output = r"""
 oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
                                  ---  PhyML 20110223  ---                                             
                            http://www.atgc-montpellier.fr/phyml                                          
                         Copyright CNRS - Universite Montpellier II                                 
 oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

. Sequence filename: 			part1.phy
. Data set: 				#1
. Tree topology: 			fixed
. Initial tree: 			user tree (optimised_start_tree.phy)
. Model of nucleotides substitution: 	HKY85
. Number of taxa: 			16
. Log-likelihood: 			-14771.02516
. Unconstrained likelihood: 		-5085.09326
. Parsimony: 				2017
. Tree size: 				0.09873
. Discrete gamma model: 		No
  - Number of categories: 		8
  - Relative rate in class 1: 		1.00000 [prop=0.125000] 		
  - Relative rate in class 2: 		1.00000 [prop=0.125000] 		
  - Relative rate in class 3: 		1.00000 [prop=0.125000] 		
  - Relative rate in class 4: 		1.00000 [prop=0.125000] 		
  - Relative rate in class 5: 		1.00000 [prop=0.125000] 		
  - Relative rate in class 6: 		1.00000 [prop=0.125000] 		
  - Relative rate in class 7: 		1.00000 [prop=0.125000] 		
  - Relative rate in class 8: 		1.00000 [prop=0.125000] 		
. Transition/transversion ratio: 	1.649
. Nucleotides frequencies:
  - f(A)= 0.19025
  - f(C)= 0.32044
  - f(G)= 0.29719
  - f(T)= 0.19212

. Run ID:				none
. Random seed:				1299409504
. Subtree patterns aliasing:		no
. Version:				20110223
. Time used:				0h0m2s (2 seconds)

 oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
 Suggested citations:
 S. Guindon, JF. Dufayard, V. Lefort, M. Anisimova, W. Hordijk, O. Gascuel
 "New algorithms and methods to estimate maximum-likelihood phylogenies: assessing the performance of PhyML 3.0."
 Systematic Biology. 2010. 59(3):307-321.

 S. Guindon & O. Gascuel
 "A simple, fast, and accurate algorithm to estimate large phylogenies by maximum likelihood"
 Systematic Biology. 2003. 52(5):696-704.
 oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
"""

if __name__ == '__main__':
    p = Parser()
    print p.parse(typical_output)
