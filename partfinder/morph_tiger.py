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

from alignment import Alignment
import numpy as np
from config import the_config
from util import PartitionFinderError

# Create a set partition for each column in the alignment
def create_set_parts(alignment):
    log.debug("Creating set partitions")
    morph_align = alignment.data.T
    set_parts = []
    part_set_dict = {}
    for col in morph_align:
        part_set_dict = {}
        for tax,i in enumerate(col):
            if i in part_set_dict:
                part_set_dict[i].append(tax)
            elif i != ord('?') and i != ord('-'):
                part_set_dict[i] = [tax]
        interim = []
        for i in part_set_dict:
            interim.append(part_set_dict[i])
        set_parts.append(interim)
    return set_parts

# Calculate similarity between two set partitions
def axpi(set_part_1, set_part_2):
    total = len(set_part_2)
    count = 0
    for i in set_part_2:
        sub_part = False
        for j in set_part_1:
            if set(i).issubset(j):
                sub_part = True
                count += 1
                break
    return float(count)/total

# Estimate rates by comparing each set partition axpi score
def calculate_rates(set_parts):
    log.debug("Estimating TIGER rates")
    rates = []
    total = len(set_parts)
    for count0,i in enumerate(set_parts):
        number = 0
        for count1,j in enumerate(set_parts):
            if count0 == count1:
                pass
            else:
                number += axpi(i, j)
        rates.append([number/(total-1)])
    return rates


if __name__ == "__main__":
    set_parts = [[[1,3],[2],[4]], [[1],[2],[3],[4]], [[1,2,3],[4]], [[1,2],[3,4]], [[1,2,3,4]]]
    print calculate_rates(set_parts)
