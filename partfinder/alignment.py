"""Loading, Saving, Parsing and Analysing Alignment Files

    Fasta is currently the only thing supported
"""
import logging
log = logging.getLogger("alignment")

import tempfile
import os

from pyparsing import (
    Word, OneOrMore, alphas, Suppress, Optional, Group, stringEnd,
    delimitedList, ParseException, line, lineno, col, LineStart, restOfLine,
    LineEnd, White, Literal)

import config, modelgen

class AlignmentError(Exception):
    pass

class FastaParser(object):
    """Parses a fasta definition and returns species codon pairs"""
    def __init__(self):
        self.sequences = {}
        self.seqlen = None

        # Some syntax that we need, but don't bother looking at
        GREATER = Literal(">")

        # I think this covers it...
        CODONS = Word(alphas + "!*-") 

        seqname = Suppress(LineStart() + GREATER) + restOfLine
        seqname.setParseAction(lambda toks: "".join(toks))
        seqcodons = OneOrMore(CODONS + Suppress(LineEnd()))
        seqcodons.setParseAction(lambda toks: "".join(toks))

        # Any sequence is the name follow by the codons
        seq = Group(seqname("species") + seqcodons("codons"))

        # Main parser: one or more definitions 
        self.root_parser = OneOrMore(seq) + stringEnd

    def parse(self, s):
        try:
            defs = self.root_parser.parseString(s)
        except ParseException, p:
            log.error("Error in Fasta Parsing:" + str(p))
            raise AlignmentError

        # Don't make a dictionary yet, as we want to check it for repeats
        return [(x.species, x.codons) for x in defs]

# Stateless singleton parser
the_parser = FastaParser()

class Alignment(object):
    def __init__(self, name):
        self.name = name
        log.debug("Created %s", self)

    def __str__(self):
        return "Alignment(%s)" % self.name

    @property
    def source_path(self):
        return self.path + ".fasta"

    @property
    def analysis_path(self):
        return self.path + ".out"

    @property
    def path(self):
        # Defer to something that get's defined differently in subclasses
        # And cache it.
        if not hasattr(self, "_path"):
            self._path = self.get_path()
            log.debug("making path %s", self._path)
        return self._path

    def get_path(self):
        # Don't use this class -- use one of the subclasses below
        raise NotImplemented

    def from_parser_output(self, defs):
        """A series of species / sequences tuples
        e.g def = ("dog", "GATC"), ("cat", "GATT")
        """
        species = {}
        slen = None
        for spec, seq in defs: 
            log.debug("Found Sequence for %s: %s...", spec, seq[:20])
            if spec in species:
                log.error("Repeated species name '%s' is repeated "
                          "in alignment", spec)
                raise AlignmentError 

            # Assign it
            species[spec] = seq

            if slen is None:
                slen = len(seq)
            else:
                if len(seq) != slen:
                    log.error("Sequence length of %s "
                              "differs from previous sequences", spec)
                    raise AlignmentError

        self.species = species
        self.sequence_length = slen

    def read_source(self):
        if not self.exists():
            log.error("Cannot find sequence file '%s'", self.source_path)
            raise AlignmentError

        log.debug("Roading fasta file '%s'", self.source_path)
        text = open(self.source_path, 'r').read()
        self.from_parser_output(the_parser.parse(text))
        
    def write_source(self):
        fd = open(self.source_path, 'w')
        log.debug("Writing %s to fasta file '%s'", self, self.source_path)
        for species, sequence in self.species.iteritems():
            fd.write(">%s\n" % species)
            fd.write("%s\n" % sequence)

    def source_exists(self):
        return os.path.exists(self.source_path)
    def analysis_exists(self):
        return os.path.exists(self.analysis_path)

    def check_against_saved(self):
        # TODO: check that it is the same...?
        pass

    def analyse(self):
        if self.source_exists():
            # We're already written a file of this name
            log.debug("%s already found at '%s'", self, self.source_path)
            self.check_against_saved()
        else:
            # Otherwise write it out
            self.write_source()

        fresh_analysis = True
        if self.analysis_exists():
            # The analysis has already been written!
            log.debug("Reading in previous analysis of %s", self)
            output = file(self.analysis_path, 'r').read()
            fresh_analysis = False
        else:
            output = modelgen.run(self.source_path)
            log.debug("Saving analysis output of %s to %s", self,
                      self.analysis_path)
            open(self.analysis_path, 'w').write(output)

        log.debug("Parsing ModelGenerator output for %s", self)
        result = modelgen.parse(output)

        if fresh_analysis:
            # Show them how long it took.
            log.debug("New analysis of %s took %d seconds", self, result.processing_time)
        return result

class SourceAlignment(Alignment):
    """The source alignment that is found in the config folder"""
    def __init__(self, name):
        Alignment.__init__(self, name)
        self.read_source()

    def get_path(self):
        return os.path.join(config.settings.base_path, self.name)

class SubsetAlignment(Alignment):
    """Created an alignment based on some others and a subset definition"""
    def __init__(self, name, source, subset):
        """create an alignment for this subset"""
        Alignment.__init__(self, name)
        # First, check to see if it exists already. If so, don't bother
        # creating it again.
        # self.species = species
        # self.sequence_length = slen
        # TODO Check sequence len against subset 
        species = {}

        # Pull out the columns we need
        for sname, old_seq in source.species.iteritems():
            new_seq = ''.join([old_seq[i] for i in subset.columns])
            species[sname] = new_seq
        self.species = species

    def get_path(self):
        return os.path.join(config.settings.output_path, self.name)

class TestAlignment(Alignment):
    """Good for testing stuff"""
    def __init__(self, name, text):
        Alignment.__init__(self, name)
        self.from_parser_output(the_parser.parse(text))

    def get_path(self):
        return os.path.join(config.settings.output_path, self.name)

if __name__ == '__main__':
    test_alignment = r"""
>spp1
CTTGAGGTTCAGAATGGTAATGAA------GTGCTGG
>spp2
CTTGAGGTACAAAATGGTAATGAG------AGCCTGG
>spp3
CTTGAGGTACAGAATAACAGCGAG------AAGCTGG
>spp4
CTCGAGGTGAAAAATGGTGATGCT------CGTCTGG
    """
    logging.basicConfig(level=logging.DEBUG)
    config.initialise("~/tmp")
    a = TestAlignment('test', test_alignment)
    res = a.analyse()
    # print res.AIC



