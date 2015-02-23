# Copyright (C) 2012 Robert Lanfear and Brett Calcott
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details. You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# PartitionFinder also includes the PhyML program, the RAxML program, and the
# PyParsing library, all of which are protected by their own licenses and
# conditions, using PartitionFinder implies that you agree with those licences
# and conditions as well.

"""Loading, Saving, Parsing Alignment Files

    See the phyml details here:
    http://www.atgc-montpellier.fr/phyml/usersguide.php?type=command

"""
import logtools
log = logtools.get_logger()

import os
from pyparsing import (Word, OneOrMore, alphas, nums, Suppress, stringEnd,
                       ParseException, restOfLine, LineEnd, ZeroOrMore, Upcase)
from util import PartitionFinderError
import numpy
import array

class AlignmentError(PartitionFinderError):
    pass


class AlignmentParser(object):
    """Parses an alignment and returns species sequence tuples"""

    # I think this covers it...
    BASES = Upcase(Word(alphas + "?.-"))

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
        seq_start = sequence_name("species") + bases(
            "sequence") + Suppress(LineEnd())
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
            self.root_parser.parseString(s)

        except ParseException as p:
            log.error("Error in Alignment Parsing:" + str(p))
            log.error("A common cause of this error is having whitespace"
                      ", i.e. spaces or tabs, in the species names. Please check this and remove"
                      " all whitespace from species names, or replace them with e.g. underscores")

            raise AlignmentError

        # Check that all the sequences are equal length
        slen = None
        names = set()
        for nm, seq in self.sequences:
            if nm in names:
                log.error("Repeated species name '%s' is repeated "
                          "in alignment", nm)
                raise AlignmentError

            names.add(nm)
            if slen is None:
                # Use the first as the test case
                slen = len(seq)
            else:
                if len(seq) != slen:
                    log.error(
                        "Bad alignment file: Not all species have the same sequences length")
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


class Alignment(object):
    def __init__(self):
        self.species = []
        self.sequence_length = 0
        self.data = None

    @property
    def species_count(self):
        return len(self.species)

    def __str__(self):
        return "Alignment(%s species, %s codons)"\
               % (self.species_count, self.sequence_length)

    def same_as(self, other):
        if self.sequence_length != other.sequence_length:
            log.warning("Alignments not the same, length differs %s: %s",
                        self.sequence_length, other.sequence_length)
            return False

        if self.species_count != other.species_count:
            log.warning("Alignments not the same. "
                        "This alignment has %s species, the alignment from the previous "
                        "analysis had %s.", len(self.species), len(other.species))
            return False

        return True

    def parse(self, text):
        """Parse the sequence, then transfer data from the parser
        Note: parser returns tuples like ("dog", "GATC"), ("cat", "GATT")
        """
        p = AlignmentParser()
        p.parse(text)

        # Allocate a numpy array using unsigned char
        d = numpy.zeros((p.species_count, p.sequence_length), 'u1')
        self.sequence_length = p.sequence_length
        for i, (spec, codons) in enumerate(p.sequences):
            self.species.append(spec)
            # TODO: Is there a better way to do this directly into numpy?
            # Typecode "B" makes a unsigned int
            d[i] = array.array("B", codons)
        self.data = d

    def read(self, pth):
        if not os.path.exists(pth):
            log.error("Cannot find sequence file '%s'", pth)
            raise AlignmentError

        log.debug("Reading alignment file '%s'", pth)
        text = open(pth, 'rU').read()
        self.parse(text)

    def write(self, pth):
        fd = open(pth, 'w')
        log.debug("Writing phylip file '%s'", pth)
        self.write_phylip(fd)
        fd.close()

    def write_phylip(self, stream):
        species_count = len(self.species)
        stream.write("%d %d\n" % (species_count, self.sequence_length))
        for i in range(species_count):
            spec = self.species[i]
            sequence = self.data[i]
            # We use a version of phylip which can have longer species names,
            # up to 100
            shortened = "%s    " % (spec[:99])
            stream.write(shortened)
            stream.write(sequence.tostring())
            stream.write("\n")


class SubsetAlignment(Alignment):

    """Create an alignment based on some others and a subset definition"""
    def __init__(self, source, subset):
        """create an alignment for this subset"""
        Alignment.__init__(self)

        # Let's do a basic check to make sure that the specified sites
        # aren't > alignment length
        site_max = max(subset.columns) + 1
        log.debug("Max site in data_blocks: %d; max site in alignment: %d"
                  % (site_max, source.sequence_length))
        if site_max > source.sequence_length:
            log.error("Site %d is specified in [data_blocks], "
                      "but the alignment only has %d sites. "
                      "Please check." % (site_max, source.sequence_length))
            raise AlignmentError

        self.species = source.species
        # Pull out the columns we need using the magic of numpy indexing
        self.data = source.data[:,subset.columns]
        self.sequence_length = len(subset.columns)
        assert self.sequence_length == self.data.shape[1]


class TestAlignment(Alignment):

    """Good for testing stuff"""
    def __init__(self, text):
        Alignment.__init__(self)
        self.parse(text)
