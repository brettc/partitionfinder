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
import model_loader as mo
from util import memoize
from config import the_config

from model_utils import get_num_params
log = logtools.get_logger()



@memoize
def get_model_commandline(modelstring):
    '''
    Input a model string, and get the piece of the raxml command line that defines that model
    '''

    commandline = the_config.available_models.query("name=='%s'" % modelstring).raxml_commandline.values[0]
    return commandline


def get_model_difficulty(modelstring):
    """
    Input a model string like HKY+I+G or LG+G+F, and a guess about how long it
    takes to analyse Right now, this is done with a simple hack. I just return
    a number that is the number of params plus a modifier for extra stuff like
    +I and +G the hardest models are +I+G, then +G, then +I this is just used
    to rank models for ordering the analysis. The return is a 'difficulty'
    score that can be used to rank models for queing and do the hardest first
    """

    elements = modelstring.split("+")

    model_params = get_num_params(modelstring)

    difficulty = 0
    if "G" in elements[1:]:
        difficulty += 2000
    if "I" in elements[1:]:
        difficulty += 1000
    if "F" in elements[1:]:
        difficulty += 3000
    if "X" in elements[1:]:
        difficulty += 4000

    if the_config.datatype == "protein" and "GTR" in modelstring:
        # that's a tough model with 189 free parameters
        difficulty += 10000

    if "LG4" in modelstring:
        # these models are hard
        difficulty += 9000

    extras = modelstring.count("+")
    total = model_params + extras + difficulty
    log.debug("Model: %s Difficulty: %d" % (modelstring, total))

    return total

if __name__ == "__main__":
    print "  ",
    print "Name".ljust(15),
    print "Params".ljust(10),
    print "Diff".ljust(10),
    print "CommandLine"
    for i, model in enumerate(get_all_dna_models()):
        print str(i+1).rjust(2), 
        print model.ljust(15),
        print str(get_num_params(model)).ljust(10),
        print str(get_model_difficulty(model)).ljust(10),
        print get_model_commandline(model)
    for i, model in enumerate(get_all_protein_models()):
        print str(i+1).rjust(2), 
        print model.ljust(15),
        print str(get_num_params(model)).ljust(10),
        print str(get_model_difficulty(model)).ljust(10),
        print get_model_commandline(model)

