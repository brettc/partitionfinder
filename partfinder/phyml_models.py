import logging
log = logging.getLogger("phyml")

import config

# TODO need some error checking!

# number of free parameters in substitution model, listed as "model+base_frequencies"
# and the model string for PhyML as the second of the tuple
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
def get_all_models():
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
def get_num_params(modelstring):
    '''
    Input a model string like HKY+I+G, and get the number of parameters
    '''
    elements = modelstring.split("+")
    model_name = elements[0]
    model_params = _base_models[model_name][0]
    extras = modelstring.count("+")
    total = model_params+extras
    return total
    
@memoize
def get_model_commandline(modelstring):
    '''
    Input a model string, and get the PhyML command line
    '''

    # This is always the same
    commandline = ["-o lr "]

    elements = modelstring.split("+")
    model_name = elements[0]
    commandline.append(_base_models[model_name][1])

    # Everything but the first element
    extras = elements[1:]

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
