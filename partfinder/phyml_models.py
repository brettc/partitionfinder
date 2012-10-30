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
# and the model string for PhyML as the second of the tuple.
_base_models = {
    "JC"    :   (0+0, "-m 000000 -f '0.25, 0.25, 0.25, 0.25'"),
    "K80"   :   (1+0, "-m 010010 -f '0.25, 0.25, 0.25, 0.25'"),
    "TrNef" :   (2+0, "-m 010020 -f '0.25, 0.25, 0.25, 0.25'"),
    "K81"   :   (2+0, "-m 012210 -f '0.25, 0.25, 0.25, 0.25'"),
    "TVMef" :   (4+0, "-m 012314 -f '0.25, 0.25, 0.25, 0.25'"),
    "TIMef" :   (3+0, "-m 012230 -f '0.25, 0.25, 0.25, 0.25'"),
    "SYM"   :   (5+0, "-m 012345 -f '0.25, 0.25, 0.25, 0.25'"),
    "F81"   :   (0+3, "-m 000000 -f e"),
    "HKY"   :   (1+3, "-m 010010 -f e"),
    "TrN"   :   (2+3, "-m 010020 -f e"),  
    "K81uf" :   (2+3, "-m 012210 -f e"),
    "TVM"   :   (4+3, "-m 012314 -f e"),
    "TIM"   :   (3+3, "-m 012230 -f e"),
    "GTR"   :   (5+3, "-m 012345 -f e")
}

# number of free parameters in substitution model, listed as "aa_frequencies"
# and the model string for PhyML as the second of the tuple
_base_protein_models = {
    "LG"            :   (0, "-m LG        -d aa"),
    "WAG"           :   (0, "-m WAG       -d aa"),
    "mtREV"         :   (0, "-m mtREV     -d aa"),
    "Dayhoff"       :   (0, "-m Dayhoff   -d aa"),
    "DCMut"         :   (0, "-m DCMut     -d aa"),
    "JTT"           :   (0, "-m JTT       -d aa"),
    "VT"            :   (0, "-m VT        -d aa"),
    "Blosum62"      :   (0, "-m Blosum62  -d aa"),
    "CpREV"         :   (0, "-m CpREV     -d aa"),
    "RtREV"         :   (0, "-m RtREV     -d aa"),
    "MtMam"         :   (0, "-m MtMam     -d aa"),
    "MtArt"         :   (0, "-m MtArt     -d aa"),
    "HIVb"          :   (0, "-m HIVb      -d aa"),
    "HIVw"          :   (0, "-m HIVw      -d aa"),
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
def get_all_dna_models():
    '''
    Return a list of all implemented _base_models
    '''
    model_list = []
    for model in _base_models.keys():
        model_list.append(model)
        model_list.append("%s+I"   %(model))
        model_list.append("%s+G"   %(model))
        model_list.append("%s+I+G" %(model))
    return model_list

@memoize
def get_all_protein_models():
    '''
    Return a list of all implemented _base__protein_models
    '''
    model_list = []
    for model in _base_protein_models.keys():
        model_list.append(model)
        model_list.append("%s+F"     %(model))
        model_list.append("%s+I"     %(model))
        model_list.append("%s+G"     %(model))
        model_list.append("%s+I+G"   %(model))
        model_list.append("%s+I+F"   %(model))
        model_list.append("%s+G+F"   %(model))
        model_list.append("%s+I+G+F" %(model))
    return model_list

@memoize
def get_mrbayes_models():
    '''
    Return a list of all models implemented in MrBayes. Thanks to Ainsley Seago for this.
    '''
    mrbayes_base_models = ["JC", "F81", "K80", "HKY", "SYM", "GTR"]
    model_list = []
    for model in mrbayes_base_models:
        model_list.append(model)
        model_list.append("%s+I"   %(model))
        model_list.append("%s+G"   %(model))
        model_list.append("%s+I+G" %(model))
    return model_list

def get_beast_models():
    '''
    Return a list of all models implemented in BEAST v1.7.2.
    '''
    beast_base_models = ["K80", "TrNef", "SYM", "HKY", "TrN", "GTR"]
    model_list = []
    for model in beast_base_models:
        model_list.append(model)
        model_list.append("%s+I"   %(model))
        model_list.append("%s+G"   %(model))
        model_list.append("%s+I+G" %(model))
    return model_list


@memoize
def get_raxml_models():
    '''
    Return a list of all models implemented in RaxML. Thanks to Ainsley Seago for this.
    '''
    model_list = ["GTR+G", "GTR+I+G"]
    return model_list

@memoize
def get_protein_models():
    '''
    Return a list of all protein models implemented in PhyML
    '''
    model_list = [
		"LG",
		"cheese"
	]
    return model_list



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
 
 
 
@memoize
def get_model_commandline(modelstring):
    '''
    Input a model string, and get the PhyML command line
    '''

    # This is always the same - optimise brlens and model, not tree
    commandline = ["-o lr "]

    elements = modelstring.split("+")
    model_name = elements[0]

    # Everything but the first element
    extras = elements[1:]

    if model_name in _base_models.keys(): #DNA models
        commandline.append(_base_models[model_name][1])
    else: #protein models
        commandline.append(_base_protein_models[model_name][1])
        if "F" in extras:
            commandline.append("-f e") #emprical AA frequencies (+19 params)
        else:
            commandline.append("-f m") #AA frequences from the model (+0 params)


    if "I" in extras:
        commandline.append("-v e")
    if "G" in extras:
        commandline.append("-a e")
        commandline.append("-c 4")
    else:
        commandline.append("-c 1")

    return " ".join(commandline)

if __name__ == "__main__":
    print "  ",
    print "Name".ljust(12),
    print "Params".ljust(10),
    print "CommandLine"
    for i, model in enumerate(get_all_models()):
        print str(i+1).rjust(2), 
        print model.ljust(12),
        print str(get_num_params(model)).ljust(10),
        print get_model_commandline(model)
    for model in get_protein_models():
        print model

