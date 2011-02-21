"""Loading and Saving Alignment Files

    Fasta is currently the only thing supported
"""
import logging
log = logging.getLogger("alignment")

from pyparsing import (
    Word, OneOrMore, alphas, Suppress, Optional, Group, stringEnd,
    delimitedList, ParseException, line, lineno, col, LineStart, restOfLine,
    LineEnd, White, Literal)

import modelgen
import tempfile

class AlignmentError(Exception):
    pass

class FastaParser(object):
    """Loads the Parser files and validates them"""
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
            log.error(p.format_message())
            raise AlignmentError

        # Don't make a dictionary, as we want to check it for repeats
        return [(x.species, x.codons) for x in defs]

the_parser = FastaParser()

class Alignment(object):
    def __init__(self, *defs):
        """A series of species / sequences

        e.g 
        
        Alignment(("dog", "GATC"), ("cat", "GATT"))
        """

    def from_parser_output(self, defs):
        species = {}
        sequence_len = None
        for spec, seq in defs: 
            log.debug("Found Sequence for %s: %s...", spec, seq[:20])
            if spec in species:
                log.error("Repeated species name '%s' is repeated "
                          "in alignment", spec)
                raise AlignmentError 

            if sequence_len is None:
                sequence_len = len(seq)
            else:
                if len(seq) != sequence_len:
                    log.error("Sequence length of %s "
                              "differs from previous sequences", spec)
                    raise AlignmentError

        self.species = species
        self.sequence_length = sequence_len

    def from_file(self):
        pass

    def write(self, fd):
        # f = open(pth, 'w')
        # log.debug("Writing fasta file '%s'", pth)
        for species, sequence in self.species.iteritems():
            fd.write(">%s\n" % species)
            fd.write("%s\n" % sequence)

    def exists(self):
        pass

    def from_columns(self, source, columns):
        pass

    def path(self):
        raise NotImplemented

    def analyse(self):
        

        pass
        # modelgen stuff here

class SubsetAlignment(Alignment):
    """Created on the fly, if it doesn't already exist"""
    def __init__(self, unique_name, columns):
        # First, check to see if it exists already. If so, don't bother
        # creating it again.
        #

        pass

    def path(self):
        return config.settings.output_path


class SourceAlignment(Alignment):
    """Read from file"""
    def __init__(self, source_name):
        pass

    def path(self):
        return config.settings.base_path

class TestAlignment(SourceAlignment):
    """Good for testing stuff"""
    def __init__(self, text):
        self.name = "unknown"
        self.saved = False
        defs = the_parser.parse(text)
        self.from_parser_output(defs)

    @property
    def path(self):
        if not self.saved:
            self._tempfile = tempfile.NamedTemporaryFile()
            log.debug("Writing alignment to %s", self._tempfile.name)
            self.write(self._tempfile.file)
            self.saved = True

        return self._tempfile.name


if __name__ == '__main__':
    import sys
    import modelgen
    import config
    logging.basicConfig(level=logging.DEBUG)
    config.initialise("~/tmp")
    a = TestAlignment(r"""
>spp1
CTTGAGGTTCAGAATGGTAATGAA------GTGCTGG
>spp2
CTTGAGGTACAAAATGGTAATGAG------AGCCTGG
>spp3
CTTGAGGTACAGAATAACAGCGAG------AAGCTGG
>spp4
CTCGAGGTGAAAAATGGTGATGCT------CGTCTGG
""")

    modelgen.run(a.path)






