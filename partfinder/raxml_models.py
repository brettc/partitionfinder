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

import logging
log = logging.getLogger("analysis")

import config

# TODO need some error checking!

# number of free parameters in substitution model, listed as "model+base_frequencies"
_base_models = {
    "GTR"   :   (5+3, "")
}

# number of free parameters in substitution model, listed as "aa_frequencies"
_base_protein_models = {
    "DAYHOFF"   :   (0, ""),
    "DCMUT"     :   (0, ""),
    "JTT"       :   (0, ""),
    "MTREV"     :   (0, ""),
    "WAG"       :   (0, ""),
    "RTREV"     :   (0, ""),
    "CPREV"     :   (0, ""),
    "VT"        :   (0, ""),
    "BLOSUM62"  :   (0, ""),
    "MTMAM"     :   (0, ""),
    "LG"        :   (0, ""),
 }

# All the functions in here return the same thing with the same parameters, 
# this just caches the return ...
def memoize(f):
    cache= {}
    def memf(*x):
        if x not in cache:
            cache[x] = f(*x)
        return cache[x]
    return memf

@memoize
def get_protein_models_gamma():
    '''
    Return a list of all implemented _base__protein_models in RAxML
    NB there are NO models in RAxML without Gamma
    '''
    model_list = []
    for model in _base_protein_models.keys():
        model_list.append("%s+G"     %(model))
        model_list.append("%s+G+F"     %(model))
    return model_list

@memoize
def get_protein_models_gammaI():
    '''
    Return a list of all implemented _base__protein_models in RAxML with invariant sites
    '''
    model_list = []
    for model in _base_protein_models.keys():
        model_list.append("%s+I+G"     %(model))
        model_list.append("%s+I+G+F"    %(model))
    return model_list

def get_all_protein_models():
    model_list = get_protein_models_gamma() + get_protein_models_gammaI()
    return model_list

@memoize
def get_dna_models_gamma():
    '''
    Just one model in RAxML with +G.
    '''
    model_list = ["GTR+G"]
    return model_list

@memoize
def get_dna_models_gammaI():
    '''
    Just one model in RAxML with I+G.
    '''
    model_list = ["GTR+I+G"]
    return model_list

@memoize
def get_all_dna_models():
    model_list = get_dna_models_gamma() + get_dna_models_gammaI()
    return model_list

@memoize
def get_all_models():
    model_list = get_all_DNA_models() + get_all_protein_models()
    return model_list

@memoize
def get_model_commandline(modelstring):
    '''
    Input a model string, and get the piece of the raxml command line that defines that model
    '''
    commandline = '-m '
    elements = modelstring.split("+")
    model_name = elements[0]

    # Everything but the first element
    extras = elements[1:]

    if model_name in _base_models.keys(): #DNA models
        commandline = ''.join([commandline, "GTRGAMMA"])
        if "I" in extras:
            commandline = ''.join([commandline, "I"])
    else: #protein models, look like this 'PROTGAMMAILGF
        commandline = ''.join([commandline, "PROTGAMMA"])
        if "I" in extras:
            commandline = ''.join([commandline, "I"])
        commandline = ''.join([commandline, model_name])        
        if "F" in extras:
            commandline = ''.join([commandline, "F"])

    return commandline


@memoize
def get_num_params(modelstring):
    '''
    Input a model string like HKY+I+G or LG+G+F, and get the number of parameters
    '''
    elements = modelstring.split("+")
    model_name = elements[0]
    if model_name in _base_models.keys():
        model_params = _base_models[model_name][0]
    else:
        model_params = _base_protein_models[model_name][0]
        if "F" in elements[1:]:
            model_params = model_params+19-1 #the -1 here is to account for the fact we add 1 for the + in '+F' below
    
    extras = modelstring.count("+")
    total = model_params+extras
    log.debug("Model: %s Params: %d" %(modelstring, total))

    return total

@memoize
def get_model_difficulty(modelstring):
    '''
    Input a model string like HKY+I+G or LG+G+F, and a guess about how long it takes to analyse
    Right now, this is done with a simple hack. I just return a number that is the number of params
    plus a modifier for extra stuff like +I and +G
    the hardest models are +I+G, then +G, then +I
    this is just used to rank models for ordering the analysis
    The return is a 'difficulty' score that can be used to rank models
    '''
    elements = modelstring.split("+")

    model_params = get_num_params(modelstring)
    
    difficulty = 0
    if "G" in elements[1:]:
        difficulty = difficulty + 2000
    if "I" in elements[1:]: 
        difficulty = difficulty + 1000
    
    extras = modelstring.count("+")
    total = model_params+extras+difficulty
    log.debug("Model: %s Difficulty: %d" %(modelstring, total))

    return total

if __name__ == "__main__":
    print "  ",
    print "Name".ljust(15),
    print "Params".ljust(10),
    print "Diff".ljust(10),
    print "CommandLine"
    for i, model in enumerate(get_all_DNA_models()):
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

