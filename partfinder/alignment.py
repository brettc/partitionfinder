#Copyright (C) 2011 Robert Lanfear and Brett Calcott
#
#This program is free software: you can redistribute it and/or modify it
#under the terms of the GNU General Public License as published by the
#Free Software Foundation, either version 3 of the License, or (at your
#option) any later version.
#
#This program is distributed in the hope that it will be useful, but
#WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#General Public License for more details. You should have received a copy
#of the GNU General Public License along with this program.  If not, see
#<http://www.gnu.org/licenses/>. PartitionFinder also includes the PhyML
#program and the PyParsing library both of which are protected by their
#own licenses and conditions, using PartitionFinder implies that you
#agree with those licences and conditions as well.

"""Loading, Saving, Parsing Alignment Files

    See the phyml details here:
    http://www.atgc-montpellier.fr/phyml/usersguide.php?type=command

"""
import logging
log = logging.getLogger("alignment")

import os

from pyparsing import (
    Word, OneOrMore, alphas, nums, Suppress, Optional, Group, stringEnd,
    delimitedList, ParseException, line, lineno, col, LineStart, restOfLine,
    LineEnd, White, Literal, Combine, Or, MatchFirst, ZeroOrMore)

from util import PartitionFinderError
class AlignmentError(PartitionFinderError):
    pass

class AlignmentParser(object):
    """Parses an alignment and returns species sequence tuples"""
    
    # I think this covers it...
    BASES = Word(alphas + "?.-")

    def __init__(self):
        self.sequence_length = None
        self.species_count = None
        self.sequences = []
        self.current_sequence = 0

        self.root_parser = self.phylip_parser() + stringEnd

    def phylip_parser(self):

        INTEGER = Word(nums) 
        INTEGER.setParseAction(lambda x: int(x[0]))

        header = INTEGER("species_count") + INTEGER("sequence_length") +\
                Suppress(restOfLine)
        header.setParseAction(self.set_header)

        sequence_name = Word(
            alphas + nums + "!#$%&\'*+-./;<=>?@[\\]^_`{|}~", 
            max=100)

        # Take a copy and disallow line breaks in the bases
        bases = self.BASES.copy()
        bases.setWhitespaceChars(" \t")
        seq_start = sequence_name("species") + bases("sequence") + Suppress(LineEnd())
        seq_start.setParseAction(self.set_seq_start)
        seq_start_block = OneOrMore(seq_start)
        seq_start_block.setParseAction(self.set_start_block)

        seq_continue = bases("sequence") + Suppress(LineEnd())
        seq_continue.setParseAction(self.set_seq_continue)

        seq_continue_block = Suppress(LineEnd()) + OneOrMore(seq_continue)
        seq_continue_block.setParseAction(self.set_continue_block)

        return header + seq_start_block + ZeroOrMore(seq_continue_block)

    def set_header(self, text, loc, tokens):
        self.sequence_length = tokens.sequence_length
        self.species_count = tokens.species_count

    def set_seq_start(self, text, loc, tokens):
        self.sequences.append([tokens.species, tokens.sequence])
        self.current_sequence += 1

    def set_start_block(self, tokens):
        # End of block
        # Reset the counter
        self.current_sequence = 0

    def set_seq_continue(self, text, loc, tokens):
        append_to = self.sequences[self.current_sequence]
        append_to[1] += tokens.sequence
        self.current_sequence += 1

    def set_continue_block(self, tokens):
        self.current_sequence = 0

    def parse(self, s):
        try:
            defs = self.root_parser.parseString(s)
        except ParseException, p:
            log.error("Error in Alignment Parsing:" + str(p))
            log.error("A common cause of this error is having whitespace"
            ", i.e. spaces or tabs, in the species names. Please check this and remove"
            " all whitespace from species names, or replace them with e.g. underscores")
                        
            raise AlignmentError

        # Check that all the sequences are equal length
        slen = None
        for nm, seq in self.sequences:
            if slen is None:
                # Use the first as the test case
                slen = len(seq)
            else:
                if len(seq) != slen:
                    log.error("Bad alignment file: Not all species have the same sequences length")
                    raise AlignmentError

        # Not all formats have a heading, but if we have one do some checking
        if self.sequence_length is not None:
            if self.sequence_length != slen:
                log.error("Bad Alignment file: sequence length count in header does not match"
                " sequence length in file, please check")
                raise AlignmentError

        if self.species_count is not None:
            if len(self.sequences) != self.species_count:
                log.error("Bad Alignment file: species count in header does not match" 
                " number of sequences in file, please check")
                raise AlignmentError

        return self.sequences

def parse(s):
    return AlignmentParser().parse(s)

class Alignment(object):
    def __init__(self):
        self.species = {}
        self.sequence_len = 0

    def __str__(self):
        return "Alignment(%s species, %s codons)" % self.species, self.sequence_len

    def same_as(self, other):
        return self.sequence_len == other.sequence_len and self.species == other.species

    def from_parser_output(self, defs):
        """A series of species / sequences tuples
        e.g def = ("dog", "GATC"), ("cat", "GATT")
        """
        species = {}
        sequence_len = None
        for spec, seq in defs: 
            # log.debug("Found Sequence for %s: %s...", spec, seq[:20])
            if spec in species:
                log.error("Repeated species name '%s' is repeated "
                          "in alignment", spec)
                raise AlignmentError 

            # Assign it
            species[spec] = seq

            if sequence_len is None:
                sequence_len = len(seq)
            else:
                if len(seq) != sequence_len:
                    log.error("Sequence length of %s "
                              "differs from previous sequences", spec)
                    raise AlignmentError
        log.debug("Found %d species with sequence length %d", 
                  len(species), sequence_len)

        # Overwrite these
        self.species = species
        self.sequence_len = sequence_len

    def read(self, pth):
        if not os.path.exists(pth):
            log.error("Cannot find sequence file '%s'", pth)
            raise AlignmentError

        log.info("Reading alignment file '%s'", pth)
        text = open(pth, 'rU').read()
        self.from_parser_output(parse(text))

    def write(self, pth):
        self.write_phylip(pth)

    def write_phylip(self, pth):
        fd = open(pth, 'w')
        log.debug("Writing phylip file '%s'", pth)

        species_count = len(self.species)
        sequence_len = len(iter(self.species.itervalues()).next())

        fd.write("%d %d\n" % (species_count, sequence_len))
        for species, sequence in self.species.iteritems():
            # we use a version of phylip which can have longer species names, up to 100
            shortened = "%s    " %(species[:99])
            fd.write(shortened)
            fd.write(sequence)
            fd.write("\n")
        fd.close()

class SubsetAlignment(Alignment):
    """Create an alignment based on some others and a subset definition"""
    def __init__(self, source, subset):
        """create an alignment for this subset"""
        Alignment.__init__(self)

        #let's do a basic check to make sure that the specified sites aren't > alignment length
        site_max = max(subset.columns)
        log.debug("Max site in data_blocks: %d; max site in alignment: %d" %(site_max, source.sequence_len))
        if site_max>source.sequence_len:
            log.error("Site %d is specified in [data_blocks], but the alignment only has %d sites. Please check." %(site_max, source.sequence_len)) 
            raise AlignmentError

        # Pull out the columns we need
        for species_name, old_sequence in source.species.iteritems():
            new_sequence = ''.join([old_sequence[i] for i in subset.columns])
            self.species[species_name] = new_sequence

        if not self.species:
            log.error("No species found in %s", self)
            raise AlignmentError

        self.sequence_len = len(self.species.itervalues().next())

class TestAlignment(Alignment):
    """Good for testing stuff"""
    def __init__(self, text):
        Alignment.__init__(self)
        self.from_parser_output(the_parser.parse(text))

