import logging
log = logging.getLogger("subset")
from fasta import write_fasta
import os

# from os import path, mkdir
# from MinimalSubsets import process_partition_file
# from string import strip, join
# from alignIO import import_fasta_as_dict
# from run_modelselection import run_modelgenerator, extract_results_modelgenerator
#

class Subset(object):
    """A set of alignment columns 

    These are the building blocks of partitioning schemes - each one has a
    given set of alignment columns as well as it's own alignment file, and
    output files from e.g. modelgenerator, each subset has likelihood scores,
    numbers of parameters, etc. etc.  some subsets will be made up of a single
    minimal subset, but most will be made up of >1 minimal subset
    """
    def __init__(self, partitions, subset_id):
        self.partitions = partitions
        self.subset_id = subset_id
        
        # Append all of the columns in the partition
        self.columns = []
        self.columnset = set()
        for part_name in subset_id:
            part = self.partitions[part_name]
            self.columns += part.columns
            self.columnset |= part.columnset
        self.columns.sort()

    def process(self, config):
        """Do all the processing"""
        self.write_alignment(config)

    @property
    def fname(self):
        s = ['[%s]' % sid for sid in sorted(list(self.subset_id))]
        return ''.join(s) + ".fasta"

        # self.allminimalsubsets    = allminimalsubsets #a dict of all of the MinimalSubset objects, which identify the partitions
        # self.input_aln_path   = input_aln_path    #the filepath of the main alignment that this subset is derived from
        # self.description      = description       #a description of the subset supplied by the user, as a string
        # self.identifier           = None              #a set: the description turned into a set of the minimal subset names that it contains
        # self.minimalsubsets   = None              #a list of minimal subsets that comprise this Subset
        # self.alignment            = None              #the alignment filepath for this subset
        # self.outputfile           = None              #the outputfile from e.g. modlegenerator for this subset
        # self.lnL              = None              #a dictionary of lnL results from different models
        # self.AIC              = None              #a dictionary of AIC results from different models
        # self.AICc             = None              #a dictionary of AICc results from different models
        # self.AIC              = None              #a dictionary of BIC results from different models
        # self.models               = None              #a list of the models we wish to consider for this alignment - usually this will be all possible models, but sometimes, e.g. if we're using MrBayes after the analysis, we only want to consider those models which MrBayes can implement
        # self.lnL_best         = None              #this subset's model name for best log likelihood score
        # self.AIC_best         = None              #this subset's AIC best model (as long as the best model is in the list of models)
        # self.AICc_best            = None              #this subset's AICc best model (as long as the best model is in the list of models)
        # self.BIC_best         = None              #this subset's BIC best model (as long as the best model is in the list of models)
        # self.parse_description()
        # self.get_minimalsubsets()
        # self.create_alignment()
        # self.analyse_alignment()
        # self.parse_results()
        # self.get_bests()
    
    def write_alignment(self, config):
        """create an alignment for this subset"""
        align = {}

        align_path = os.path.join(config.output_path, self.fname)
        if os.path.exists(align_path):
            log.debug(
                "Fasta file '%s' already exists, not rewriting",
                os.path.basename(align_path))
            return 

        # Pull out the columns we need
        for species, old_seq in config.sequence.iteritems():
            new_seq = ''.join([old_seq[i] for i in self.columns])
            align[species] = new_seq

        write_fasta(align_path, align)
        self.align_path = align_path

    # def analyse_alignment(self):
        # """analyse the alignment for a subset"""
        # thisfile_path = path.dirname(path.abspath(__file__))
        # modelgenerator_path = "%s/../programs/modelgenerator.jar" %(thisfile_path)
        # 
        # output_filename = "%s_modelgenerator.out" %(self.alignment)

        # print output_filename

        # if not(path.isfile(output_filename)):
            # self.outputfile = run_modelgenerator(modelgenerator_path, self.alignment)
        # else:
            # self.outputfile = output_filename
            
    # def parse_results(self):
        # """parse_results from a modelgenerator output file"""
        # self.AIC, self.AICc, self.BIC, self.lnL = extract_results_modelgenerator(self.outputfile)
        # 
    # def get_bests(self):
        # """get best scores according to the input list of models"""
        # #if model list wasn't specified, use all models output by whatever the analysis program was
        # if self.models == None:
            # self.models = self.lnL.keys()
        # #pick initial values for scores
        # first_model = self.models[0]
        # AIC_best = self.AIC[first_model]+1
        # AICc_best = self.AICc[first_model]+1
        # BIC_best = self.BIC[first_model]+1
        # lnL_best = self.lnL[first_model]-1
        # for model in self.models:
            # if self.AIC[model]<AIC_best: 
                # self.AIC_best = model
                # AIC_best = self.AIC[model]
            # if self.AICc[model]<AICc_best: 
                # self.AICc_best = model
                # AICc_best = self.AICc[model]
            # if self.BIC[model]<BIC_best: 
                # self.BIC_best = model
                # BIC_best = self.BIC[model]
            # if self.lnL[model]>lnL_best: 
                # self.lnL_best = model
                # lnL_best = self.lnL[model]
        # 
            # 
