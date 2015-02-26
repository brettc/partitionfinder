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

import logtools
import pandas as pd
import os
from config import the_config

log = logtools.get_logger()
from util import PartitionFinderError

def get_num_params(modelstring):
    """
    Input a model string like HKY+I+G or LG+G+F, and get the number of
    parameters
    """

    m = the_config.available_models.query("name=='%s'" %modelstring).matrix_params.values[0]
    b = the_config.available_models.query("name=='%s'" %modelstring).basefreq_params.values[0]
    r = the_config.available_models.query("name=='%s'" %modelstring).ratevar_params.values[0]

    total = m+b+r

    log.debug("Model: %s Params: %d" % (modelstring, total))

    return total



