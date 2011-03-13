"""Loading, Saving, Parsing Alignment Files
"""
import logging
log = logging.getLogger("alignment")

import tempfile
import os

from pyparsing import (
    Word, OneOrMore, alphas, nums, Suppress, Optional, Group, stringEnd,
    delimitedList, ParseException, line, lineno, col, LineStart, restOfLine,
    LineEnd, White, Literal, Combine, Or, MatchFirst)

import config

# No longer using fasta, but we'll keep it around for the moment...
# Should really detect which it is...
# alignment_format = 'fasta'
alignment_format = 'phy'

class AlignmentError(Exception):
    pass

class AlignmentParser(object):
    """Parses a fasta definition and returns species codon pairs"""
    
    # I think this covers it...
    CODONS = Word(alphas + "!*-") 

    def __init__(self):
        self.sequences = {}
        self.seqlen = None

        # TODO: We should be able to read both of these...!
        if alignment_format == 'phy':
            self.root_parser = self.phylip_parser() + stringEnd
        elif alignment_format == 'fasta':
            self.root_parser = self.fasta_parser() + stringEnd

    def fasta_parser(self):
        # Some syntax that we need, but don't bother looking at
        GREATER = Literal(">")

        seqname = Suppress(LineStart() + GREATER) + restOfLine
        seqname.setParseAction(lambda toks: "".join(toks))
        seqcodons = OneOrMore(self.CODONS + Suppress(LineEnd()))
        seqcodons.setParseAction(lambda toks: "".join(toks))

        # Any sequence is the name follow by the codons
        seq = Group(seqname("species") + seqcodons("codons"))
        sequences = OneOrMore(seq)
        # Main parser: one or more definitions 
        return sequences("sequences")

    def phylip_parser(self):

        INTEGER = Word(nums) 
        INTEGER.setParseAction(lambda x: int(x[0]))

        header = Group(INTEGER("species_count") +
                       INTEGER("codon_count") + Suppress(restOfLine))

        seqname = Word(alphas + nums, max=10)

        # Take a copy and disallow line breaks in the codons
        codons = self.CODONS.copy()
        codons.setWhitespaceChars(" \t")
        codonrepeats = OneOrMore(codons)

        codonrepeats.setParseAction(lambda x: ''.join(x))
        seq = Group(seqname("species") + codonrepeats("codons")) + Suppress(LineEnd())
        sequences = OneOrMore(seq)
        return header("header") + sequences("sequences")

    def parse(self, s):
        try:
            defs = self.root_parser.parseString(s)
        except ParseException, p:
            log.error("Error in Alignment Parsing:" + str(p))
            raise AlignmentError
        print defs

        # if we have a header, do some checking
        if defs.header:
            if len(defs.sequences) != defs.header.species_count:
                log.error("Bad Alignment file: species count does not match")
                raise AlignmentError

            if len(defs.sequences[0][1]) != defs.header.codon_count:
                log.error("Bad Alignment file: codon count does not match")
                raise AlignmentError

        # Don't make a dictionary yet, as we want to check it for repeats
        return [(x.species, x.codons) for x in defs.sequences]

# Stateless singleton parser
the_parser = AlignmentParser()

class Alignment(object):
    def __init__(self, name):
        self.name = name
        log.debug("Created %s", self)

    def __str__(self):
        return "Alignment(%s)" % self.name

    def from_parser_output(self, defs):
        """A series of species / sequences tuples
        e.g def = ("dog", "GATC"), ("cat", "GATT")
        """
        species = {}
        slen = None
        for spec, seq in defs: 
            # log.debug("Found Sequence for %s: %s...", spec, seq[:20])
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
        log.debug("Found %d species with sequence length %d", 
                  len(species), slen)
        return species, slen

    def read(self, pth):
        if not os.path.exists(pth):
            log.error("Cannot find sequence file '%s'", pth)
            raise AlignmentError

        log.debug("Reading alignment file '%s'", pth)
        text = open(pth, 'r').read()
        return self.from_parser_output(the_parser.parse(text))

    def write(self, pth):
        if alignment_format == 'phy':
            self.write_phylip(pth)
        elif alignment_format is 'fasta':
            self.write_fasta(pth)

        log.error("Undefined Alignment Format")
        raise AlignmentError

    def write_fasta(self, pth):
        fd = open(pth, 'w')
        log.debug("Writing %s to fasta file '%s'", self, pth)
        for species, sequence in self.species.iteritems():
            fd.write(">%s\n" % species)
            fd.write("%s\n" % sequence)

    def write_phylip(self, pth):
        fd = open(pth, 'w')
        log.debug("Writing %s to phylip file '%s'", self, pth)

        species_count = len(self.species)
        codon_count = len(iter(self.species.itervalues()).next())

        fd.write("%d %d\n" % (species_count, codon_count))
        for species, sequence in self.species.iteritems():
            # Species is max 10 long, clip it and fill with spaces
            shortened = species[:9].ljust(10, ' ')
            fd.write(shortened)
            fd.write(sequence)
            fd.write("\n")

class SourceAlignment(Alignment):
    """The source alignment that is found in the config folder"""
    def __init__(self, name):
        Alignment.__init__(self, name)
        # self.read_source()

    # def get_path(self):
        # return os.path.join(config.settings.base_path, self.name)

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
        s, l = self.from_parser_output(the_parser.parse(text))
        self.species = s
        self.seqlen = l

    # def get_path(self):
        # return os.path.join(config.settings.output_path, self.name)

def test_phylip():
    logging.basicConfig(level=logging.DEBUG)
    config.initialise("~/tmp")
    test_alignment = r"""16 13
A GCCGGCGGGGAAA
	B	GCCGGCGGGGAAA
	C	GCCGGCGGGGAAG
	D	GCCAGCGGGGAAG
E	GCCCGGGGGGAAG
	F	GCCCGTGGGGAAG
	G	GCCGGCGAGGAAG
	H	GCCAGCGGGGAAA
	I	GCCGGCGGGGAAG
	J	GCC GGCGGGGAAG
	K	GCCGGCGGGGAAT
	L	GCCGGCGGGGAAA
	M	GTCGGCGGGGAAA
	N	GCCGGCGGGGAAA
	O	GCCGGCGGGGAAA
	P	GCC GGCGGGGAAT
"""
    print the_parser.parse(test_alignment)

def test_fasta():
    test_alignment = r"""
>spp1
ATTGAGGTTCAGAATGGTAATGAA------GTGCTGG
>spp2
CTTGAGGTACAAAATGGTAATGAG------AGCCTGG
>spp3
CTTGAGGTACAGAATAACAGCGAG------AAGCTGG
>spp4
ATCGAGGTGAAAAATGGTGATGCT------CGTCTGG
    """
    logging.basicConfig(level=logging.DEBUG)
    config.initialise("~/tmp")
    a = TestAlignment('test', test_alignment)
    print a
    # res = a.analyse()
    # print res.AIC

if __name__ == '__main__':
    test_phylip()
