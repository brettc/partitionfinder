import logging
log = logging.getLogger("modelgen")

import config

import subprocess, shlex, tempfile, os, shutil

from pyparsing import (
    Word, OneOrMore, alphas, nums, Suppress, Group, stringEnd, ParseException,
    line, lineno, col, LineStart, SkipTo, LineEnd,
    )

class ModelGeneratorError(Exception):
    pass

# TODO: get num gamma cat from default settings
def run(alignment_path, num_gamma_cat=4):

    # Run it in a temporary directory, so we can ditch all the output. We read
    # the output directly from the program via stdout, so we don't need any of
    # the files it writes, they are just extraneous cruft
    output_dir = tempfile.mkdtemp()
    log.debug("Created temporary folder %s for modelgenerator output", output_dir)
    mg = config.settings.modelgen_path
    command = "java -jar %s %s %d" % (mg, alignment_path, num_gamma_cat)
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
    finally:
        log.debug("Removing temporary folder %s", output_dir)
        shutil.rmtree(output_dir)

    if p.returncode != 0:
        log.error("Modelgenerator program failed to execute successfully")
        raise ModelGeneratorError

    # Otherwise this should contain the output. It will need parsing
    return stdout

class ModelGeneratorResult(object):
    def __init__(self):
        self.AIC = {}
        self.AICc = {}
        self.BIC = {}
        self.lnL = {}

    def __str__(self):
        return "ModelGeneratorResult()"

# Make the parser a singleton
class Parser(object):
    def __init__(self):

        # A line of repeating -----
        HYPHENS = LineStart() + Word("-") + LineEnd()

        INTEGER = Word(nums).setParseAction(lambda x: int(x[0]))
        FLOAT = Word(nums + '.-').setParseAction(lambda x: float(x[0]))
        MODEL = Word(alphas + nums + '+')

        # Skip everything till we get to the table (2 lines of hyphens)
        # Not sure why "include = False" does not work in SkipTo...
        preamble = Suppress(
            SkipTo(HYPHENS) + HYPHENS + 
            SkipTo(HYPHENS) + HYPHENS
        )

        result_group = Group(FLOAT("result") + MODEL("model") + FLOAT("lnl"))
        table_row = Group(LineStart() + Suppress(INTEGER) + 
                     result_group("AIC") + 
                     result_group("AICc") +
                     result_group("BIC") + 
                     Suppress(LineEnd()))

        # Final timing info...
        timing = Group(Suppress("Analysis took:") + INTEGER("min") +
                       Suppress("minutes") + INTEGER("sec") +
                       Suppress("seconds"))
        table = Group(OneOrMore(table_row))

        self.root_parser = \
                preamble + table("results") + timing("timing") + stringEnd

    def parse(self, input_string):
        try:
            parsed = self.root_parser.parseString(input_string)
        except ParseException, p:
            log.error(p.format_message())
            raise ModelGeneratorError

        # Fill out one of the results
        r = ModelGeneratorResult()
        for line in parsed.results:
            r.AIC[line.AIC.model] = line.AIC.result
            r.AICc[line.AICc.model] = line.AICc.result
            r.BIC[line.BIC.model] = line.BIC.result
            r.lnL[line.AIC.model] = line.AIC.lnl

        r.processing_time = parsed.timing.min * 60 + parsed.timing.sec
        return r
    
# The parser can be a singleton (as it is stateless)
the_parser = Parser()

def parse(input_string):
    log.debug("Parsing ModelGenerator output....")
    return the_parser.parse(input_string)

if __name__ == '__main__':
    import tempfile, os, shutil
    logging.basicConfig(level=logging.DEBUG)
    align = r""">spp1
CTTGAGGTTCAGAATGGTAATGAA------GTGCTGG
>spp2
CTTGAGGTACAAAATGGTAATGAG------AGCCTGG
>spp3
CTTGAGGTACAGAATAACAGCGAG------AAGCTGG
>spp4
CTCGAGGTGAAAAATGGTGATGCT------CGTCTGG
"""
    config.initialise("~/tmp")
    outname = os.path.join(config.settings.output_path, 'temp.fasta')
    f = open(outname, 'wb')
    f.write(align)
    f.close()
    outp = run(outname)
    # outp = open('output.txt', 'r').read()
    result = parse(outp)
    print result.AIC
    print result.processing_time


