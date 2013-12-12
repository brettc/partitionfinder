# Copyright (C) 2012-2013 Robert Lanfear and Brett Calcott
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

import logtools
log = logtools.get_logger(__file__)

import os
import numpy
import tables

import raxml_models
import phyml_models

int_type = numpy.int32
float_type = numpy.float32


def _model_string_maxlen():
    """Calculate the field size needed for model ids"""
    all_models = \
        raxml_models.get_all_dna_models() +\
        raxml_models.get_all_protein_models() +\
        phyml_models.get_all_dna_models() +\
        phyml_models.get_all_protein_models()

    lengths = [len(m) for m in all_models]
    return max(lengths)

# TODO: This is overkill, we never need the enums. It should just be a
# dictionary
def enum(sequential):
    """
    Modified from here
    http://stackoverflow.com/questions/36932/how-can-i-represent-an-enum-in-python
    We use these values to index into the database fields
    """
    n = len(sequential)
    enums = dict(zip(sequential, range(n)))
    # We add an extra one
    enums['MAX'] = n

    # Good for looking up by variable
    enums['get'] = staticmethod(lambda x: enums[x])

    return type('Enum', (), enums)


# TODO: This should probably be in Raxml.py....
# The number of codons
codon_types = list("ATCG")
Freqs = enum(codon_types)

# Define the number of rates we record
# TODO: This should be generated using itertools.combinations(x, 2)
rates_types = ["{}_{}".format(f, t) for f, t in zip('AAACCG', 'CGTGTT')]
Rates = enum(rates_types)


def make_results_datatype():
    # 32 is the md5 length
    subset_id_length = 32
    model_id_length = _model_string_maxlen()

    layout = [
        ('subset_id', 'S{}'.format(subset_id_length)),
        ('model_id', 'S{}'.format(model_id_length)),
        ('seconds', int_type),
        ('params', int_type),
    ]

    # Now add the floating point fields
    flds = "lnl alpha aic aicc bic site_rate".split()
    for f in flds:
        layout.append((f, float_type))

    # Now add frequencies and rates
    # Note that these are added as embedded in an extra dimension
    # We can access them via the Enums above.
    # EG. Given a numpy record X, we can use X['freqs'][Freqs.A]
    layout.extend([
        ('freqs', float_type, (Freqs.MAX)),
        ('rates', float_type, (Rates.MAX)),
    ])

    # Now construct the numpy datatype that gives us the layout
    return numpy.dtype(layout)


# Calculate this now, it is fixed, and not analysis dependent
results_dtype = make_results_datatype()


class Database(object):

    def __init__(self, cfg):
        self.cfg = cfg
        self.path = os.path.join(self.cfg.subsets_path, 'data.db')
        if os.path.exists(self.path):
            self.h5 = tables.open_file(self.path, 'a')
            self.results = self.h5.root.results
            assert isinstance(self.results, tables.Table)
            assert self.results.indexed
        else:
            # Compression is good -- and faster, according to the pytables docs...
            f = tables.Filters(complib='blosc', complevel=5)
            self.h5 = tables.open_file(self.path, 'w', filters=f)
            self.results = self.h5.create_table(
                '/', 'results', results_dtype)
            self.results.cols.subset_id.create_csindex()

    def get_results_for_subset(self, subset):
        conditions = {'current_id':  subset.name}
        matching = self.results.read_where(
            'subset_id == current_id', conditions)
        return matching

    def save_result(self, subset, n):
        # We have to take a slice here, as pytables can't handle single
        # elements
        self.results.append(subset.result_array[n:n+1])
        self.cfg.database.results.flush()

    def close(self):
        self.h5.close()
