"""Simple Fasta Parser
"""
import logging
log = logging.getLogger("fasta")

from pyparsing import (
    Word, OneOrMore, alphas, Suppress, Optional, Group, stringEnd,
    delimitedList, ParseException, line, lineno, col, LineStart, restOfLine,
    LineEnd, White)

class FastaError(Exception):
    pass

class FastaParser(object):
    """Loads the Parser files and validates them"""
    def __init__(self):
        self.make_syntax()
        self.sequences = {}
        self.seqlen = None

    def make_syntax(self):

        # Some syntax that we need, but don't bother looking at
        GREATER = Suppress(">")

        # I think this covers it...
        CODONS = Word(alphas + "!*-") 

        seqname = LineStart() + GREATER + restOfLine
        seqcodons = OneOrMore(CODONS + Suppress(LineEnd()))
        # A parse action to collect all the lines
        seqcodons.setParseAction(self.join_codons_action)

        # Any sequence is the name follow by the codons
        seq = seqname + seqcodons
        seq.setParseAction(self.set_sequence_action)

        # Main parser: one or more definitions 
        self.fasta = OneOrMore(seq) + stringEnd

    def join_codons_action(self, text, loc, tokens):
        # Simply join them together and send it back
        return "".join(tokens)

    def set_sequence_action(self, text, loc, tokens):
        seqname = tokens[0]
        sequence = tokens[1]
        log.debug("Found Sequence for %s: %s...", seqname, sequence[:20])

        if seqname in self.sequences:
            log.error("Repeated species name '%s' at line %d", seqname,
                      lineno(loc, text))
            raise FastaError

        if self.seqlen is None:
            self.seqlen = len(sequence)
        else:
            if len(sequence) != self.seqlen:
                log.error("Sequence length of %s at line %d "
                          "differs from previous sequences", seqname,
                      lineno(loc, text))
                raise FastaError

        # We'll create it as we go
        self.sequences[seqname] = sequence

    def parse_file(self, fname):
        s = open(fname, 'r').read()
        return self.parse_configuration(s)

    def parse_configuration(self, s):
        self.fasta.parseString(s)
        return self.sequences

def read_fasta(fname):
    p = FastaParser()
    log.debug("Parsing fasta file '%s'", fname)
    return p.parse_file(fname)

def write_fasta(fname, fasta_dict):
    f = open(fname, 'w')
    log.debug("Writing fasta file '%s'", fname)
    for species, sequence in fasta_dict.iteritems():
        f.write(">%s\n" % species)
        f.write("%s\n" % sequence)
        

if __name__ == '__main__':
    test1 = r"""
>spp1
CTTGAGGTTCAGAATGGTAATGAA------GTGCTGG
>spp2
CTTGAGGTACAAAATGGTAATGAG------AGCCTGG
>spp3
CTTGAGGTACAGAATAACAGCGAG------AAGCTGG
>spp4
CTCGAGGTGAAAAATGGTGATGCT------CGTCTGG
"""
    import sys
    logging.basicConfig(level=logging.DEBUG)
    c = FastaParser()
    c.parse_configuration(test1)







