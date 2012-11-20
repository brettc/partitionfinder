from cluster import minkowski_distance, genmatrix

def sum_matrices(mat1, mat2, mat3):
    """sum three matrices of the same size"""
    result = []
    rowindex=0
    for row in mat1:
        result.append([sum(x) for x in zip(mat1[rowindex], mat2[rowindex], mat3[rowindex])] )
        rowindex += 1
    return result

def normalise_and_weight(matrix, normaliser, weight):
    """normalise and weight a 2D matrix (list of lists) using the passed value"""
    
    if normaliser==0:
        #then our matrix must have had a maximum distance of zero anyway, indicating that 
        #everything is already identical (e.g. same amino acid model)
        #in this case we don't need to do any work, because these distances won't count
        #when we sum the matrices anyway, so just return without worrying
        return matrix
        
    rowindex=0
    for row in matrix:
        colindex = 0 # keep track of where we are in the matrix
        for cell in row:
            matrix[rowindex][colindex] = cell*float(weight)/float(normaliser)
            colindex += 1
        rowindex += 1

    return matrix
    
    
def get_minmax(matrix):
    """return the min and max distances in a 2D list
    """

    maxdistance=None    
    mindistance=None
    rowindex=0
    for row in matrix:
        colindex = 0 # keep track of where we are in the matrix
        for cell in row:
            #ignore anything on or below the diagonal
            if (colindex > rowindex):
                if (cell<mindistance or mindistance is None): mindistance=cell
                if (cell>maxdistance or maxdistance is None): maxdistance=cell
            colindex += 1
        rowindex += 1

    return mindistance, maxdistance

def get_closest_pair(matrix, subsets):
    """return the closest pair of subsets defined by a matrix 
    """

    mindistance=None
    rowindex=0
    for row in matrix:
        colindex = 0 # keep track of where we are in the matrix
        for cell in row:
            #ignore anything on or below the diagonal
            if (colindex > rowindex):
                if (cell<mindistance or mindistance is None): 
                    mindistance=cell
                    min_row = rowindex
                    min_col = colindex
            colindex += 1
        rowindex += 1

    sub1 = subsets[min_row]
    sub2 = subsets[min_col]

    return sub1, sub2


def get_closest_subsets(scheme, weights):
    """Find the two closest subsets in a scheme
    """
    #1. get the parameter lists for each subset
    subsets = [] #a list of subset names, so we know the order things appear in the list
    rates = [] #tree length
    freqs = [] #amino acid or base frequencies
    model = [] #model parameters e.g. A<->C
    
    for s in scheme.subsets:
        param_dict = s.get_param_values()
        subsets.append(s)
        rates.append([param_dict["rate"]])
        freqs.append(param_dict["freqs"])
        model.append(param_dict["model"])

    #2. for each parameter we get a matrix of euclidean distances
    rates_matrix = genmatrix(list=rates, combinfunc=minkowski_distance, symmetric=True, diagonal=0)
    freqs_matrix = genmatrix(list=freqs, combinfunc=minkowski_distance, symmetric=True, diagonal=0)
    model_matrix = genmatrix(list=model, combinfunc=minkowski_distance, symmetric=True, diagonal=0)


    #3. Normalise and weight those euclidean distances,
    min, max = get_minmax(rates_matrix)
    rates_matrix = normalise_and_weight(rates_matrix, max, weights["rate"])
    min, max = get_minmax(freqs_matrix)
    freqs_matrix = normalise_and_weight(freqs_matrix, max, weights["freqs"])
    min, max = get_minmax(model_matrix)
    model_matrix = normalise_and_weight(model_matrix, max, weights["model"])
    
    #4. sum the matrices
    distance_matrix = sum_matrices(rates_matrix, freqs_matrix, model_matrix)

    #5. get the closest pair
    sub1, sub2 = get_closest_pair(distance_matrix, subsets)
    
    return sub1, sub2

def get_nearest_neighbour_scheme(
        start_scheme, scheme_name, cfg, weights = {"rate": 1, "freqs": 1, "model": 1}):
    """The idea here is to take a scheme, and perform some analyses to find a neighbouring
    scheme, where the neighbour has one less subset than the current scheme. 
    Really this is just progressive clustering, but specified to work well with PartitionFinder
    
    The weights argument allows us to assign different weights to different model parameters
    
    """
        
    import subset
    import scheme
    
    #1. First we get the closest subsets, based on some weights
    #   the weights are [tree_size, model_params, base_freqs]
    sub1, sub2 = get_closest_subsets(start_scheme, weights)

    #2. Next we create a new subset that merges those two subsets
    newsub_parts = list(sub1.partitions) + list(sub2.partitions)
    newsub = subset.Subset(*tuple(newsub_parts))
    
    #3. Then we define a new scheme with those merged subsets
    subs = [s for s in start_scheme.subsets]

    #pop out the two subsets we're going to join together
    subs.remove(sub1)
    subs.remove(sub2)

    subs.append(newsub)

    scheme = (scheme.Scheme(cfg, str(scheme_name), *tuple(subs)))

    return scheme    
        
    