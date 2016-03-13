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
log = logtools.get_logger()

import os
import numpy
import tables
from itertools import combinations

import raxml_models
import phyml_models

int_type = numpy.int32
float_type = numpy.float32


def _model_string_maxlen():
    """Calculate the field size needed for model ids"""
    # hardcoded for convenience. Could be dynamically set in future.
    # the current longest is: BLOSUM62+I+G+X, i.e. 14 chars. 
    # so we just over double it, for safety

    return 30


class DataLayout(object):
    def __init__(self, letters=None):
        self.letters = letters
        if letters is not None:
            self.make_results_and_freqs()
        else:
            # We just fake an entry, 
            self.letter_indexes = { 'EMPTY': 0 }
            self.rate_indexes = { 'EMPTY': 0 }
            self.letter_size = 1
            self.rate_size = 1

        self.data_type = self.make_datatype()

    def make_results_and_freqs(self):
        l = list(self.letters) 
        self.letter_indexes = dict(zip(l, range(len(l))))
        self.letter_size = len(self.letter_indexes)

        ri = {}
        for i, rate in enumerate(combinations(l, 2)):
            # We need both directions as either could be used to look it up
            # ie. A <-> C or C <-> A
            f, t = rate
            ri["%s_%s" % (f, t)] = i
            ri["%s_%s" % (t, f)] = i

        self.rate_indexes = ri
        self.rate_size = len(ri) / 2 

    def get_empty_record(self):
        return numpy.zeros(1, self.data_type)

    def make_datatype(self):
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

        # Now add frequencies and rate. These are added as embedded in an extra dimension

        layout.extend([
            ('freqs', float_type, self.letter_size),
            ('rates', float_type, self.rate_size),
        ])

        # Now construct the numpy datatype that gives us the layout
        return numpy.dtype(layout)


class DataRecord(object):
    def __init__(self, cfg):
        self.__dict__['_data'] = cfg.data_layout.get_empty_record()

    def __getattr__(self, name):
        return self._data[name]

    def __setattr__(self, name, value):
        self._data[name] = value

    def __str__(self):
        return "DataRecord<lnl:%s, tree_size:%s, secs:%s>" % (
            self.lnl, self.site_rate, self.seconds)


class Database(object):

    def __init__(self, cfg):
        self.cfg = cfg
        self.path = os.path.join(self.cfg.subsets_path, 'data.db')
        self.results = None
        if os.path.exists(self.path):
            try:
                self.h5 = tables.open_file(self.path, 'a')
                self.results = self.h5.root.results
            except:
                # If anything fails, we just create a new database...
                log.warning("""Failed to open existing database at %s, or
                database is corrupted. Creating a new one""", self.path)
                self.results = None

        # Something went wrong!
        if not self.results:
            try:
                # Try closing this, just in case
                self.h5.close()
            except:
                pass

            # Compression is good -- and faster, according to the pytables docs...
            f = tables.Filters(complib='blosc', complevel=5)
            self.h5 = tables.open_file(self.path, 'w', filters=f)
            self.results = self.h5.create_table(
                '/', 'results', cfg.data_layout.data_type)
            self.results.cols.subset_id.create_csindex()

        assert isinstance(self.results, tables.Table)
        assert self.results.indexed

    def get_results_for_subset(self, subset):
        conditions = {'current_id':  subset.subset_id}
        matching = self.results.read_where(
            'subset_id == current_id', conditions)
        return matching

    def is_empty(self):
        return self.results.nrows == 0

    def save_result(self, subset, n):
        # We have to take a slice here, as pytables can't handle single
        # elements
        self.results.append(subset.result_array[n:n+1])
        self.cfg.database.results.flush()

    def close(self):
        self.h5.close()
