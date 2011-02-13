"""Simple Fasta Parser
"""
import logging
log = logging.getLogger("fasta")

from pyparsing import (
    Word, OneOrMore, alphas, Suppress, Optional, Group, stringEnd,
    delimitedList, ParseException, line, lineno, col, LineStart, restOfLine,
    LineEnd, White)

class ParserError(Exception):
    def __init__(self, text, loc, msg):
        """Used for our own parsing problems"""
        self.line = line(loc, text)
        self.col = col(loc, text)
        self.lineno = lineno(loc, text)
        self.msg = msg

    def format_message(self):
        return "%s at line:%s, column:%s" % (self.msg, self.lineno, self.col)
        
class Parser(object):
    """Loads the Parser files and validates them"""
    def __init__(self):
        super(Parser, self).__init__()
        self.make_syntax()
        self.sequences = {}

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
        log.debug("Found Sequence for %s" % seqname)

        if seqname in self.sequences:
            raise ConfigurationError(text, loc, "Repeated sequence name")

        # We'll create it as we go
        self.sequences[seqname] = sequence

    def parse_file(self, fname):
        log.info("Reading File %s", fname)
        s = open(fname, 'r').read()
        return self.parse_configuration(s)

    def parse_configuration(self, s):
        try:
            # We ignore the return as we build everything in the actions
            self.fasta.parseString(s)
        except ConfigurationError, c:
            log.error(c.format_message())
            self.error = True

def parse_file(fname):
    p = Parser()
    return p.parse_file(fname)



def test_basic():
    test1 = r"""
>spp1
CTTGAGGTTCAGAATGGTAATGAA------GTGCTGGT
GCTGGAAGTTCAGCAGCAGCTCGGCGGCGG
>spp2
CTT

GAGGTACAAAATGGTAATGAG------AGCCTGGTG
>spp3
CTTGAGGTACAGAATAACAGCGAG------AAGCTGGT
>spp4

CTCGAGGTGAAAAATGGTGATGCT------CGTCTGGT
"""
    import sys
    logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
    c = Parser()
    c.parse_configuration(test1)
    print c.sequences

if __name__ == '__main__':
    test_basic()
    # print c.parts







