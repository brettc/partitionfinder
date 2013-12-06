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
import tables

import raxml_models, phyml_models

def _model_string_maxlen():
    """Calculate the field size needed for model ids"""
    all_models = \
            raxml_models.get_all_dna_models() +\
            raxml_models.get_all_protein_models() +\
            phyml_models.get_all_dna_models() +\
            phyml_models.get_all_protein_models()

    lengths = [len(m) for m in all_models]
    return max(lengths)

# Record description -- add and remove fields here...
class ResultDescription(tables.IsDescription):
    subset_id = tables.StringCol(32)
    model_id = tables.StringCol(_model_string_maxlen())
    lnl = tables.Float32Col()
    seconds = tables.Int32Col()
    alpha = tables.Float32Col()
    params = tables.Int32Col()
    aic = tables.Float32Col()
    bic = tables.Float32Col()
    aicc = tables.Float32Col()
    site_rate = tables.Float32Col()

# This is the numpy equivalent, we'll use it to create
results_dtype = tables.dtype_from_descr(ResultDescription)

class Database(object):
    def __init__(self, cfg):
        self.cfg = cfg
        self.path = os.path.join(self.cfg.subsets_path, 'data.db')
        if os.path.exists(self.path):
            self.load()
        else:
            self.create()

    def create(self):
        # f = tables.Filters(complib='blosc', complevel=5)
        self.h5 = tables.openFile(self.path, 'w') # , filters=f)
        self.results = self.h5.createTable(
            '/', 'results', ResultDescription)

    def load(self):
        self.h5 = tables.openFile(self.path, 'a')

    def save_result(self):
        pass
