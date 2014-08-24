import csv
import sys

def likelihood_parser(fasttree_logfile):
    list_of_dicts = []
    # Try opening the fasttree logfile
    try:
        with open(str(fasttree_logfile)) as fasttree_logfile:
            first_line = fasttree_logfile.next()
            hello = first_line.split(' ')
            # Make sure that the logfile looks like you think it should
            if hello[0] == "Command:":
                # Find the line that starts listing the likelihoods
                while True:
                    new_line = fasttree_logfile.next()
                    split_line = new_line.split('\t')
                    print split_line[0]
                    first_word = split_line[0]
                    print first_word

                    # Read in likelihoods
                    if first_word == "Gamma20LogLk":
                        list_of_dicts = list(csv.DictReader(fasttree_logfile,
                            delimiter = '\t', skipinitialspace = True))
                        break
            else:
                raise IOError('Likelihood file corrupted')
    except IOError:
        raise IOError('Could not locate likelihood file!')

    # Make a list containing all of the log likelihoods
    lnl_list = [lnl['LogLk'] for lnl in list_of_dicts]
    # Get rid of the junk at the bottom of the file that doesn't contain lnLs
    lnl_list = lnl_list[:-6]
    return lnl_list


if __name__ == '__main__':
    fasttree_logfile = sys.argv[1]
    list_of_dicts = likelihood_parser(fasttree_logfile)
    print list_of_dicts