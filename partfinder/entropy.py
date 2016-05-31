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

# defining a new function 'entropy_calc' which takes as input a 1D array p
# copied from here: http://nbviewer.ipython.org/url/atwallab.cshl.edu/teaching/QBbootcamp3.ipynb
def entropy_calc(p):
    p=p[p!=0] # modify p to include only those elements that are not equal to 0
    return np.dot(-p,np.log2(p)) # the function returns the entropy result


def get_morph_entropies(alignment):
    morph_align = alignment.data.T
    column_counts = []
    for col in morph_align:
        new_col = col[col != ord('-')]
        new_col = new_col[new_col != ord('?')]
        column_counts.append(np.unique(new_col, return_counts=True)[1])
    column_entropy = []
    for col in column_counts:
        sum = np.sum(col)
        props = np.array([num/float(sum) for num in col])
        column_entropy.append(entropy_calc(props))
    column_entropy = np.array(column_entropy)
    column_entropy = column_entropy.reshape(len(column_entropy), 1)

    return column_entropy

def sitewise_entropies(alignment):
    if the_config.datatype == 'DNA':
        log.debug("Calculating DNA entropies")
        dna_states = "ACGT"
        dna_list = [np.sum(alignment.data == ord(nuc), axis = 0) for nuc in list(dna_states)]
        states = np.array(dna_list, dtype=float)
    elif the_config.datatype == 'protein':
        log.debug("Calculating protein entropies")
        aa_states = "ARNDCQEGHILKMFPSTWYV"
        amino_list = [np.sum(alignment.data == ord(aa), axis = 0) for aa in list(aa_states)]
        states = np.array(amino_list, dtype = float)
    elif the_config.datatype == 'morphology':
        column_entropy = get_morph_entropies(alignment)
    else:
        log.error("Unknown datatype '%s'" % the_config.datatype)
        raise PartitionFinderError

    if the_config.datatype != 'morphology':
        states = states.T
        totals = np.sum(states, axis=1)
        totals.shape = len(states),1

        # for a column of all gaps, we'll have a zero total, so we just hack that here
        totals = np.where(totals==0, 1, totals)

        prob = states/totals

        column_entropy = [[entropy_calc(t)] for t in prob]
        column_entropy = np.array(column_entropy)

    return column_entropy

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def get_replacement_sites(column_entropy, columns):
    """take an array of entropies and return the nearest site
    to each site that has zero entropy"""

    zeros = np.where(column_entropy==0)[0]
    nonzeros = np.nonzero(column_entropy)[0]

    if len(nonzeros)==0:
        log.warning("The entropy of every site in your alignment is zero, cannot reassign entropies")
        return(column_entropy)
        
    # figure out the replacement values
    replacements = {} # key: zero_entropy_site; value: replacement
    for z in zeros:
        nearest_nonzero = find_nearest(nonzeros, z)
        #replacements.setdefault(columns[nearest_nonzero], []).append(columns[z])
        replacements[columns[z]] = columns[nearest_nonzero]

    return(replacements)

def sitewise_entropies_scaled(alignment): 
    """This function will calculate entropies for DNA based on the assumption
    that the states in the column are the only possible states, i.e. it
    doesn't assume four states for each column
    """
    dna_align = alignment.data.T
    column_counts = []
    for col in dna_align:
        new_col = col[col != ord('-')]
        new_col = new_col[new_col != ord('?')]
        column_counts.append(np.unique(new_col, return_counts=True)[1])
    column_entropy = []
    for col in column_counts:
        sum = np.sum(col)
        props = np.array([num/float(sum) for num in col])
        column_entropy.append([entropy_calc(props)]) 
    column_entropy = np.array(column_entropy)
    column_entropy = column_entropy.reshape(len(column_entropy), 1)
    return column_entropy



