from pyparsing import Word, Dict, OneOrMore, alphas, nums, \
        Suppress, Optional, Group, stringEnd, delimitedList, \
        pythonStyleComment, ParseException, line, lineno, col, \
        LineStart, restOfLine, LineEnd, White


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

class ConfigurationError(Exception):
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

    def make_syntax(self):

        # Some syntax that we need, but don't bother looking at
        GREATER = Suppress(">")

        # Partition Parsing
        seqname = LineStart() + GREATER + restOfLine
        seqcodons = OneOrMore(Word("AGCT-") +
                              Suppress(LineEnd())).setParseAction(self.sequence_action)

        seq = Group(seqname + seqcodons)

        self.fasta = OneOrMore(seq) + stringEnd

    def sequence_action(self, tokens):
        return "".join(tokens)

    def parse_file(self, fname):
        s = open(fname, 'r').read()
        return self.parse_configuration(s)

    def parse_configuration(self, s):
        try:
            # We ignore the return as we build everything in the actions
            return self.fasta.parseString(s)
        except ConfigurationError, c:
            print c.format_message()
        except ParseException, p:
            # print p.line
            # print p.col, p.lineno
            print str(p)

def parse_file(fname):
    p = Parser()
    return p.parse_file(fname)

if __name__ == '__main__':
    c = Parser()
    print dict(c.parse_configuration(test1).asList())
    # print c.parts







