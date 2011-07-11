"""This is where everything comes together, and we do the analysis"""

import logging
log = logging.getLogger("analysis")

import os, shutil

from alignment import Alignment, SubsetAlignment
import phyml
import threadpool
import scheme
import algorithm
import subset
import submodels

from util import PartitionFinderError
class AnalysisError(PartitionFinderError):
    pass

def make_dir(pth):
    if os.path.exists(pth):
        if not os.path.isdir(pth):
            log.error("Cannot create folder '%s'", pth)
            raise AnalysisError
    else:
        os.mkdir(pth)

class Analysis(object):
    """Performs the analysis and collects the results"""
    def __init__(self, cfg, force_restart, threads=1):
        cfg.validate()

        log.info("Beginning Analysis")
        self.threads = threads
        self.cfg = cfg
        if force_restart:
            if os.path.exists(self.cfg.output_path):
                log.warning("Deleting all previous workings in '%s'", 
                            self.cfg.output_path)
                shutil.rmtree(self.cfg.output_path)

        #check for old analyses to see if we can use the old data
        self.cfg.check_for_old_config()

        # Make some folders for the analysis
        make_dir(self.cfg.output_path)
        self.make_output_dir('subsets')
        self.make_output_dir('schemes')
        self.make_output_dir('phyml')
        self.make_output_dir('start_tree')

        self.make_alignment(cfg.alignment_path)
        self.make_tree()
        self.subsets_analysed_set = set() #a counter for user info
        self.subsets_analysed = 0 #a counter for user info
        self.total_subset_num = None
        self.schemes_analysed = 0 #a counter for user info
        self.total_scheme_num = None

    def do_analysis(self):
        if self.cfg.search == 'all':
            self.analyse_all_possible(self.cfg.models)
        elif self.cfg.search == 'user':
            self.analyse_current_schemes(self.cfg.models)
        elif self.cfg.search == 'greedy':
            self.analyse_greedy(self.cfg.models, self.cfg.model_selection)
        else:
            log.error("Search algorithm '%s' is not yet implemented", 
                        self.cf.search)
            raise AnalysisError

    def make_alignment(self, source_alignment_path):
        # Make the alignment 
        self.alignment = Alignment()
        self.alignment.read(source_alignment_path)

        # We start by copying the alignment
        self.alignment_path = os.path.join(self.cfg.output_path, 'start_tree', 'source.phy')
        if os.path.exists(self.alignment_path):
            # Make sure it is the same
            old_align = Alignment()
            old_align.read(self.alignment_path)
            if not old_align.same_as(self.alignment):
                log.error("Alignment file has changed since previous run. "
                          "You need to use the force-restart option.")
                raise AnalysisError

        else:
            self.alignment.write(self.alignment_path)

    def make_tree(self):
        # Begin by making a filtered alignment, containing ONLY those columns
        # that are defined in the subsets
        subset_with_everything = subset.Subset(*list(self.cfg.partitions))
        self.filtered_alignment = SubsetAlignment(self.alignment, 
                                                  subset_with_everything)
        self.filtered_alignment_path = os.path.join(self.cfg.output_path,
                                                    'start_tree',
                                                    'filtered_source.phy')
        self.filtered_alignment.write(self.filtered_alignment_path)

        # Now we've written this alignment, we need to lock everything in
        # place, no more adding partitions, or changing them from now on.
        self.cfg.partitions.check_against_alignment(self.alignment)
        self.cfg.partitions.finalise()

        # Now check for the tree
        tree_path = phyml.make_tree_path(self.filtered_alignment_path)
        if not os.path.exists(tree_path):
            tree_path = phyml.make_tree(self.filtered_alignment_path)
        self.tree_path = tree_path
        log.info("BioNJ tree with GTR+I+G brlens is stored here: %s", self.tree_path) 

    def make_output_dir(self, name):
        new_path = os.path.join(self.cfg.output_path, name)
        make_dir(new_path)
        setattr(self, name+"_path", new_path)

    def analyse_subset(self, sub, models):
        """Analyse the subset using the models given
        This is the core place where everything comes together
        The results are placed into subset.result
        """
        #keep people informed about what's going on
        #if we don't know the total subset number, we can usually get it like this
        if self.total_subset_num == None:
            self.total_subset_num = len(sub._cache)
        old_num_analysed = self.subsets_analysed
        self.subsets_analysed_set.add(sub.name)
        self.subsets_analysed = len(self.subsets_analysed_set)
        if self.subsets_analysed>old_num_analysed: #we've just analysed a subset we haven't seen yet
            percent_done = float(self.subsets_analysed)*100.0/float(self.total_subset_num)
            log.info("Analysing subset %d/%d: %.2f%s done" %(self.subsets_analysed,self.total_subset_num, percent_done, r"%"))

        sub_bin_path = os.path.join(self.subsets_path, sub.name + '.bin')
        # We might have already saved a bunch of results, try there first
        if not sub.results:
            sub.read_binary_summary(sub_bin_path)

        # First, see if we've already got the results loaded. Then we can
        # shortcut all the other checks
        models_done = set(sub.results.keys())
        models_required = set(models)
        models_to_do = models_required - models_done

        # Empty set means we're done
        if not models_to_do:
            log.debug("Using results that are already loaded %s", sub)
            return


        log.debug("About to analyse %s using models %s", sub, ", ".join(list(models)))

        # Make an Alignment from the source, using this subset
        sub_alignment = SubsetAlignment(self.alignment, sub)
        sub_path = os.path.join(self.subsets_path, sub.name + '.phy')
        # Add it into the sub, so we keep it around
        sub.alignment_path = sub_path

        # Maybe it is there already?
        if os.path.exists(sub_path):
            log.debug("Found existing alignment file %s", sub_path)
            old_align = Alignment()
            old_align.read(sub_path)

            # It had better be the same!
            if not old_align.same_as(sub_alignment):
                log.error("It looks like you have changed something in the "
                        "configuration and I cannot trust the old analysis. "
                        "You'll need to run the program with --force-restart")
                raise AnalysisError
        else:
            # We need to write it
            sub_alignment.write(sub_path)

        # Try and read in some previous analyses
        self.parse_results(sub, models_to_do)
        if not models_to_do:
            return

        # What is left, we actually have to analyse...
        tasks = []
        for m in models_to_do:
            a_path, out_path = phyml.make_analysis_path(self.phyml_path, sub.name, m)
            tasks.append((phyml.analyse, 
                          (m, sub_path, a_path, self.tree_path, self.cfg.branchlengths)))

        if self.threads == 1:
            self.run_models_concurrent(tasks)
        else:
            self.run_models_threaded(tasks)

        # Now parse the models we've just done
        self.parse_results(sub, models_to_do)

        # This should be empty NOW!
        if models_to_do:
            log.error("Failed to run models %s; not sure why", 
                      ", ".join(list(models_to_do)))
            raise AnalysisError

        # Now we have analysed all models for this subset, we do model selection
        sub.model_selection(self.cfg.model_selection)        
        
        # If we made it to here, we should write out the new summary
        sub_summary_path = os.path.join(self.subsets_path, sub.name + '.txt')
        sub.write_summary(sub_summary_path)
        # We also need to update this
        sub.write_binary_summary(sub_bin_path)

    def parse_results(self, sub, models_to_do):
        """Read in the results and parse them"""
        models_done = []
        for m in list(models_to_do):
            # sub.alignment_path
            a_path, out_path = phyml.make_analysis_path(self.phyml_path, sub.name, m)
            if os.path.exists(out_path):
                sub_output = open(out_path, 'rb').read()
                # Annotate with the parameters of the model
                try:
                    result = phyml.parse(sub_output)
                    sub.add_model_result(m, result)
                    # Remove the current model from remaining ones
                    models_to_do.remove(m)
                    
                    # Just used for below
                    models_done.append(m)
                except phyml.PhymlError:
                    log.warning("Failed loading parse output from %s."
                              "Output maybe corrupted. I'll run it again.",
                              out_path)

        if models_done:
            log.debug("Loaded analysis for %s, models %s", sub, ", ".join(models_done))

    def run_models_concurrent(self, tasks):
        for func, args in tasks:
            func(*args)

    def run_models_threaded(self, tasks):
        pool = threadpool.Pool(tasks, self.threads)
        pool.join()

    def analyse_scheme(self, sch, models):
        self.schemes_analysed = self.schemes_analysed + 1        
        log.info("Analysing scheme %d/%d" %(self.schemes_analysed, self.total_scheme_num))
        for sub in sch:
            self.analyse_subset(sub, models)
 
        # AIC needs the number of sequences 
        number_of_seq = len(self.alignment.species)
        sch.assemble_results(number_of_seq, self.cfg.branchlengths)
        sch.write_summary(os.path.join(self.schemes_path, sch.name+'.txt'))

    def analyse_current_schemes(self, models):
        """Process everything when search=user"""
        current_schemes = [s for s in self.cfg.schemes]
        self.total_scheme_num = len(current_schemes)
        if self.total_scheme_num>0:
            for s in current_schemes:
                 self.analyse_scheme(s, models)
            self.write_best_scheme(current_schemes)
        else:
            log.error("Search set to 'user', but no user schemes detected in .cfg file. Please check.")
            raise PartitionFinderError

    def analyse_greedy(self, models, method):
        '''A greedy algorithm for heuristic partitioning searches'''
        log.info("Performing greedy analysis")

        partnum = len(self.cfg.partitions)
        self.total_scheme_num = submodels.count_greedy_schemes(partnum)
        log.info("This will result in a maximum of %s schemes being created", self.total_scheme_num)
        self.total_subset_num = submodels.count_greedy_subsets(partnum)
        log.info("PartitionFinder will have to analyse a maximum of %d subsets of sites to complete this analysis" %(self.total_subset_num))
        if self.total_subset_num>1000000:
            log.warning("%d is a lot of subsets, this might take a long time to analyse", self.total_subset_num)
            log.warning("Perhaps consider using a different search scheme instead (see Manual)")

        #clear any schemes that are currently loaded
        # TODO Not sure we need this...
        self.cfg.schemes.clear_schemes()        
                
        #start with the most partitioned scheme
        start_description = range(len(self.cfg.partitions))
        start_scheme = scheme.create_scheme(self.cfg, 1, start_description)
        log.info("Analysing starting scheme (scheme %s)" % start_scheme.name)
        self.analyse_scheme(start_scheme, models)
        
        def get_score(my_scheme):
            #TODO: this is bad. Should use self.cfg.model_selection, or write a new method for scheme.py
            if method=="aic":
                score=my_scheme.aic
            elif method=="aicc":
                score=my_scheme.aicc
            elif method=="bic":
                score=my_scheme.bic
            else:
                log.error("Unrecognised model_selection variable '%s', please check" %(score))
                raise AnalysisError
            return score

        best_scheme = start_scheme
        best_score  = get_score(start_scheme)
                         
        round = 1
        cur_s = 2

        #now we try out all lumpings of the current scheme, to see if we can find a better one
        #and if we do, we just keep going
        while True:
            log.info("***Greedy algorithm step %d***" % round)

            #get a list of all possible lumpings of the best_scheme
            lumpings = algorithm.lumpings(start_description)

            #we reset the counters as we go, for better user information
            self.total_scheme_num = len(lumpings)
            self.schemes_analysed = 0

            best_lumping_score = None
            for lumped_description in lumpings:
                lumped_scheme = scheme.create_scheme(self.cfg, cur_s, lumped_description)
                cur_s = cur_s + 1
                self.analyse_scheme(lumped_scheme, models)
                new_score = get_score(lumped_scheme)

                if best_lumping_score==None or new_score<best_lumping_score:
                    best_lumping_score  = new_score
                    best_lumping_scheme = lumped_scheme
                    best_lumping_desc   = lumped_description
                                                
            if best_lumping_score<best_score:
                best_scheme = best_lumping_scheme
                best_score  = best_lumping_score
                start_description = best_lumping_desc				
                if len(set(best_lumping_desc))==1: #then it's the scheme with everything equal, so quit
				    break
                round = round+1

            else:
                break

        log.info("Greedy algorithm finished after %d rounds" % round)
        log.info("Highest scoring scheme is scheme %s, with %s score of %.3f" %(best_scheme.name, method, best_score))

        best_schemes_file = os.path.join(self.cfg.output_path, 'best_schemes.txt')
        best_scheme.write_summary(best_schemes_file, 'wb', "Best scheme according to Greedy algorithm, analysed with %s\n\n" % method)
        log.info("Information on best scheme is here: %s" %(best_schemes_file))

        current_schemes = [s for s in self.cfg.schemes]
        current_schemes.sort(key=lambda s: int(s.name), reverse=False)

        self.write_all_schemes(current_schemes) #this also writes a file which has info on all analysed schemes, useful for extra analysis if that's what you're interested in...

    def analyse_all_possible(self, models):

        partnum = len(self.cfg.partitions)
        self.total_scheme_num = submodels.count_all_schemes(partnum)
        log.info("Analysing all possible schemes for %d starting partitions", partnum)
        log.info("This will result in %s schemes being created", self.total_scheme_num)
        self.total_subset_num = submodels.count_all_subsets(partnum)
        log.info("PartitionFinder will have to analyse %d subsets to complete this analysis" %(self.total_subset_num))
        if self.total_subset_num>1000000:
            log.warning("%d is a lot of subsets, this might take a long time to analyse", self.total_subset_num)
            log.warning("Perhaps consider using a different search scheme instead (see Manual)")

        #clear any schemes that are currently loaded
        self.cfg.schemes.clear_schemes()

        gen_schemes = scheme.generate_all_schemes(self.cfg)
        for s in gen_schemes:
            self.analyse_scheme(s, models)
        self.write_best_scheme(gen_schemes)

    def write_best_scheme(self, list_of_schemes):
        # Which is the best?
        sorted_schemes_aic = [(s.aic, s) for s in list_of_schemes]
        sorted_schemes_aic.sort()
        sorted_schemes_aicc = [(s.aicc, s) for s in list_of_schemes]
        sorted_schemes_aicc.sort()
        sorted_schemes_bic = [(s.bic, s) for s in list_of_schemes]
        sorted_schemes_bic.sort()
        best_aic  = sorted_schemes_aic[0][1]
        best_aicc = sorted_schemes_aicc[0][1]
        best_bic  = sorted_schemes_bic[0][1]
        best_schemes_file = os.path.join(self.cfg.output_path, 'best_schemes.txt')
        best_aic.write_summary(best_schemes_file, 'wb', "Best scheme according to AIC\n")
        best_aicc.write_summary(best_schemes_file, 'ab', "\n\n\nBest scheme according to AICc\n")
        best_bic.write_summary(best_schemes_file, 'ab', "\n\n\nBest scheme according to BIC\n")
        log.info("Information on best schemes is here: %s" %(best_schemes_file))
        self.write_all_schemes(list_of_schemes) #this also writes a file which has info on all analysed schemes, useful for extra analysis if that's what you're interested in...

    def write_all_schemes(self, list_of_schemes):
        all_schemes_file = os.path.join(self.cfg.output_path, 'all_schemes.txt')
        f = open(all_schemes_file, 'wb')
        f.write("Name\tlnL\t#params\t#sites\t#subsets\tAIC\tAICc\tBIC\n")
        for s in list_of_schemes:
            f.write("%s\t%.3f\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\n" %(s.name,s.lnl,s.sum_k,s.nsites,s.nsubs,s.aic,s.aicc,s.bic))
        log.info("Information on all schemes analysed is here: %s" %(all_schemes_file))
		
    
if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    # logging.basicConfig(level=logging.INFO)
    from alignment import TestAlignment
    from partition import Partition

    # TODO: should probably reduce the size of this
    alignment = TestAlignment("""
4 2208
spp1     CTTGAGGTTCAGAATGGTAATGAA------GTGCTGGTGCTGGAAGTTCAGCAGCAGCTCGGCGGCGGTATCGTACGTACCATCGCCATGGGTTCTTCCGACGGTCTGCGTCGCGGTCTGGATGTAAAAGACCTCGAGCACCCGATCGAAGTCCCAGTTGGTAAAGCAACACTGGGTCGTATCATGAACGTACTGGGTCAGCCAGTAGACATGAAGGGCGACATCGGTGAAGAAGAGCGTTGGGCT---------------ATCCACCGTGAAGCACCATCCTATGAAGAGCTGTCAAGCTCTCAGGAACTGCTGGAAACCGGCATCAAAGTTATCGACCTGATGTGTCCGTTTGCGAAGGGCGGTAAAGTTGGTCTGTTCGGTGGTGCGGGTGTAGGTAAAACCGTAAACATGATGGAGCTTATTCGTAACATCGCGATCGAGCACTCCGGTTATTCTGTGTTTGCGGGCGTAGGTGAACGTACTCGTGAGGGTAACGACTTCTACCACGAAATGACCGACTCCAACGTTATCGAT---------------------AAAGTTTCTCTGGTTTATGGCCAGATGAACGAGCCACCAGGTAACCGTCTGCGCGTTGCGCTGACCGGTCTGACCATGGCTGAGAAGTTCCGTGACGAAGGTCGCGACGTACTGCTGTTCGTCGATAACATCTATCGTTACACCCTGGCAGGTACTGAAGTTTCAGCACTGCTGGGTCGTATGCCTTCAGCGGTAGGTTACCAGCCGACTCTGGCGGAAGAAATGGGCGTTCGCATTCCAACGCTGGAAGAGTGTGATATCTGCCACGGCAGCGGCGCTAAAGCCGGTTCGAAGCCGCAGACCTGTCCTACCTGTCACGGTGCAGGCCAGGTACAGATGCGCCAGGGCTTCTTCGCTGTACAGCAGACCTGTCCACACTGCCAGGGCCGCGGTACGCTGATCAAAGATCCGTGCAACAAATGTCACGGTCATGGTCGCGTAGAGAAAACCAAAACCCTGTCCGTAAAAATTCCGGCAGGCGTTGATACCGGCGATCGTATTCGTCTGACTGGCGAAGGTGAAGCTGGTGAGCACGGCGCACCGGCAGGCGATCTGTACGTTCAGGTGCAGGTGAAGCAGCACGCTATTTTCGAGCGTGAAGGCAACAACCTGTACTGTGAAGTGCCGATCAACTTCTCAATGGCGGCTCTTGGCGGCGAGATTGAAGTGCCGACGCTTGATGGTCGCGTGAAGCTGAAAGTTCCGGGCGAAACGCAAACTGGCAAGCTGTTCCGTATGCGTGGCAAGGGCGTGAAGTCCGTGCGCGGCGGTGCACAGGGCGACCTTCTGTGCCGCGTGGTGGTCGAGACACCGGTAGGTCTTAACGAGAAGCAGAAACAGCTGCTCAAAGATCTGCAGGAAAGTTTTGGCGGCCCAACGGGTGAAAACAACGTTGTTAACGCCCTGTCGCAGAAACTGGAATTGCTGATCCGCCGCGAAGGCAAAGTACATCAGCAAACTTATGTCCATGGTGTGCCACAGGCTCCGCTGGCGGTAACCGGTGAAACGGAAGTGACCGGTACACAGGTGCGTTTCTGGCCAAGCCACGAAACCTTCACCAACGTAATCGAATTCGAATATGAGATTCTGGCAAAACGTCTGCGCGAGCTGTCATTCCTGAACTCCGGCGTTTCCATCCGTCTGCGCGATAAGCGTGAC---GGCAAAGAAGACCATTTCCACTATGAAGGTGGTATCAAGGCGTTTATTGAGTATCTCAATAAAAATAAAACGCCTATCCACCCGAATATCTTCTACTTCTCCACCGAA---AAAGACGGTATTGGCGTAGAAGTGGCGTTGCAGTGGAACGATGGTTTCCAGGAAAACATCTACTGCTTCACCAACAACATTCCACAGCGTGATGGCGGTACTCACCTTGCAGGCTTCCGTGCGGCGATGACCCGTACGCTGAACGCTTACATGGACAAAGAAGGCTACAGCAAAAAAGCCAAA------GTCAGCGCCACCGGTGATGATGCCCGTGAAGGCCTGATTGCCGTCGTTTCCGTGAAAGTACCGGATCCGAAATTCTCCTCTCAGACTAAAGACAAACTGGTCTCTTCTGAGGTGAAAACGGCGGTAGAACAGCAGATGAATGAACTGCTGAGCGAATACCTGCTGGAAAACCCGTCTGACGCCAAAATC
spp2     CTTGAGGTACAAAATGGTAATGAG------AGCCTGGTGCTGGAAGTTCAGCAGCAGCTCGGTGGTGGTATCGTACGTGCTATCGCCATGGGTTCTTCCGACGGTCTGCGTCGTGGTCTGGAAGTTAAAGACCTTGAGCACCCGATCGAAGTCCCGGTTGGTAAAGCAACGCTGGGTCGTATCATGAACGTGCTGGGTCAGCCGATCGATATGAAAGGCGACATCGGCGAAGAAGAACGTTGGGCG---------------ATTCACCGTGCAGCACCTTCCTATGAAGAGCTCTCCAGCTCTCAGGAACTGCTGGAAACCGGCATCAAAGTTATCGACCTGATGTGTCCGTTCGCGAAGGGCGGTAAAGTCGGTCTGTTCGGTGGTGCGGGTGTTGGTAAAACCGTAAACATGATGGAGCTGATCCGTAACATCGCGATCGAACACTCCGGTTACTCCGTGTTTGCTGGTGTTGGTGAGCGTACTCGTGAGGGTAACGACTTCTACCACGAAATGACCGACTCCAACGTTCTGGAT---------------------AAAGTATCCCTGGTTTACGGCCAGATGAACGAGCCGCCGGGAAACCGTCTGCGCGTTGCACTGACCGGCCTGACCATGGCTGAGAAATTCCGTGACGAAGGTCGTGACGTTCTGCTGTTCGTCGATAACATCTATCGTTATACCCTGGCCGGTACAGAAGTATCTGCACTGCTGGGTCGTATGCCTTCTGCGGTAGGTTATCAGCCGACGCTGGCGGAAGAGATGGGCGTTCGTATCCCGACGCTGGAAGAGTGCGACGTCTGCCACGGCAGCGGCGCGAAATCTGGCAGCAAACCGCAGACCTGTCCGACCTGTCATGGTCAGGGCCAGGTGCAGATGCGTCAGGGCTTCTTCGCCGTTCAGCAGACCTGTCCGCATTGTCAGGGGCGCGGTACGCTGATTAAAGATCCGTGCAACAAATGTCACGGTCACGGTCGCGTTGAGAAAACCAAAACCCTGTCGGTCAAAATCCCGGCGGGCGTGGATACCGGCGATCGTATTCGTCTGTCAGGAGAAGGCGAAGCGGGCGAACACGGTGCACCAGCAGGCGATCTGTACGTTCAGGTCCAGGTTAAGCAGCACGCCATCTTTGAGCGTGAAGGCAATAACCTGTACTGCGAAGTGCCTATTAACTTCACCATGGCAGCCCTCGGCGGCGAGATTGAAGTCCCGACGCTGGATGGCCGGGTGAATCTCAAAGTGCCTGGCGAAACGCAAACCGGCAAACTGTTCCGCATGCGCGGTAAAGGTGTGAAATCCGTGCGCGGTGGTGCTCAGGGCGACCTGCTGTGCCGCGTGGTGGTTGAAACACCAGTCGGGCTGAACGATAAGCAGAAACAGCTGCTGAAGGACCTGCAGGAAAGTTTTGGCGGACCAACGGGCGAGAAAAACGTGGTTAACGCCCTGTCGCAGAAGCTGGAGCTGGTTATTCAGCGCGACAATAAAGTTCACCGTCAGATCTATGCGCACGGTGTGCCGCAGGCTCCGCTGGCAGTGACCGGTGAGACCGAAAAAACCGGCACCATGGTACGTTTCTGGCCAAGCTATGAAACCTTCACCAACGTTGTCGAGTTCGAATACGAGATCCTGGCAAAACGTCTGCGTGAGCTGTCGTTCCTGAACTCCGGGGTTTCTATCCGTCTGCGTGACAAGCGTGAC---GGTAAAGAAGACCATTTCCACTACGAAGGCGGCATCAAGGCGTTCGTTGAGTATCTCAATAAGAACAAAACGCCGATCCACCCGAATATCTTCTACTTCTCCACCGAA---AAAGACGGTATTGGCGTCGAAGTAGCGCTGCAGTGGAACGACGGCTTCCAGGAAAACATCTACTGCTTCACCAACAACATCCCGCAGCGCGATGGCGGTACTCACCTTGCGGGCTTCCGCGCGGCGATGACCCGTACCCTGAACGCCTATATGGACAAAGAAGGCTACAGCAAAAAAGCCAAA------GTCAGCGCTACCGGCGACGATGCGCGTGAAGGCCTGATTGCCGTTGTCTCCGTGAAGGTTCCGGATCCGAAATTCTCCTCGCAGACCAAAGACAAACTGGTCTCCTCCGAGGTGAAAACCGCGGTTGAACAGCAGATGAATGAACTGCTGAACGAATACCTGCTGGAAAATCCGTCTGACGCGAAAATC
spp3     CTTGAGGTACAGAATAACAGCGAG------AAGCTGGTGCTGGAAGTTCAGCAGCAGCTCGGCGGCGGTATCGTACGTACCATCGCAATGGGTTCTTCCGACGGTCTGCGTCGTGGTCTGGAAGTGAAAGACCTCGAGCACCCGATCGAAGTCCCGGTAGGTAAAGCGACCCTGGGTCGTATCATGAACGTGCTGGGTCAGCCAATCGATATGAAAGGCGACATCGGCGAAGAAGATCGTTGGGCG---------------ATTCACCGCGCAGCACCTTCCTATGAAGAGCTGTCCAGCTCTCAGGAACTGCTGGAAACCGGCATCAAAGTTATCGACCTGATTTGTCCGTTCGCTAAGGGCGGTAAAGTTGGTCTGTTCGGTGGTGCGGGCGTAGGTAAAACCGTAAACATGATGGAGCTGATCCGTAACATCGCGATCGAGCACTCCGGTTACTCCGTGTTTGCAGGCGTGGGTGAGCGTACTCGTGAGGGTAACGACTTCTACCACGAGATGACCGACTCCAACGTTCTGGAC---------------------AAAGTTGCACTGGTTTACGGCCAGATGAACGAGCCGCCAGGTAACCGTCTGCGCGTAGCGCTGACCGGTCTGACCATCGCGGAGAAATTCCGTGACGAAGGCCGTGACGTTCTGCTGTTCGTCGATAACATCTATCGTTATACCCTGGCCGGTACAGAAGTTTCTGCACTGCTGGGTCGTATGCCATCTGCGGTAGGTTATCAGCCTACTCTGGCAGAAGAGATGGGTGTTCGTATCCCGACGCTGGAAGAGTGTGAAGTTTGCCACGGCAGCGGCGCGAAAAAAGGTTCTTCTCCGCAGACCTGTCCAACCTGTCATGGACAGGGCCAGGTGCAGATGCGTCAGGGCTTCTTCACCGTGCAGCAAAGCTGCCCGCACTGCCAGGGCCGCGGTACCATCATTAAAGATCCGTGCACCAACTGTCACGGCCATGGCCGCGTAGAGAAAACCAAAACGCTGTCGGTAAAAATTCCGGCAGGCGTGGATACCGGCGATCGTATCCGCCTTTCTGGTGAAGGCGAAGCGGGCGAGCACGGCGCACCTTCAGGCGATCTGTACGTTCAGGTTCAGGTGAAACAGCACCCAATCTTCGAGCGTGAAGGCAATAACCTGTACTGCGAAGTGCCGATCAACTTTGCGATGGCTGCGCTGGGCGGGGAAATTGAAGTGCCGACCCTTGACGGCCGCGTTAAGCTGAAGGTACCGAGCGAAACGCAAACCGGCAAGCTGTTCCGCATGCGCGGTAAAGGCGTGAAATCCGTACGCGGTGGCGCGCAGGGCGATCTGCTGTGCCGCGTCGTCGTTGAAACTCCGGTTAGCCTGAACGAAAAGCAGAAGAAACTGCTGCGTGATTTGGAAGAGAGCTTTGGCGGCCCAACGGGGGCGAACAATGTTGTGAACGCCCTGTCCCAGAAGCTGGAGCTGCTGATTCGCCGCGAAGGCAAAACCCATCAGCAAACCTACGTGCACGGTGTGCCGCAGGCTCCGCTGGCGGTCACCGGTGAAACCGAACTGACCGGTACCCAGGTGCGTTTCTGGCCGAGCCATGAAACCTTCACCAACGTCACCGAATTCGAATATGACATCCTGGCTAAGCGCCTGCGTGAGCTGTCGTTCCTGAACTCCGGCGTCTCTATTCGCCTGAACGATAAGCGCGAC---GGCAAGCAGGATCACTTCCACTACGAAGGCGGCATCAAGGCGTTTGTTGAGTACCTCAACAAGAACAAAACCCCGATTCACCCGAACGTCTTCTATTTCAGCACTGAA---AAAGACGGCATCGGCGTGGAAGTGGCGCTGCAGTGGAACGACGGCTTCCAGGAAAATATCTACTGCTTTACCAACAACATTCCTCAGCGCGACGGCGGTACTCACCTTGCGGGCTTCCGCGCGGCGATGACCCGTACCCTGAACGCCTATATGGACAAAGAAGGCTACAGCAAAAAAGCCAAA------GTGAGCGCCACCGGTGACGATGCGCGTGAAGGCCTGATTGCCGTAGTGTCCGTGAAGGTGCCGGATCCGAAGTTCTCTTCCCAGACCAAAGACAAACTGGTTTCTTCGGAAGTGAAATCCGCGGTTGAACAGCAGATGAACGAACTGCTGGCTGAATACCTGCTGGAAAATCCGGGCGACGCAAAAATT
spp4     CTCGAGGTGAAAAATGGTGATGCT------CGTCTGGTGCTGGAAGTTCAGCAGCAGCTGGGTGGTGGCGTGGTTCGTACCATCGCCATGGGTACTTCTGACGGCCTGAAGCGCGGTCTGGAAGTTACCGACCTGAAAAAACCTATCCAGGTTCCGGTTGGTAAAGCAACCCTCGGCCGTATCATGAACGTATTGGGTGAGCCAATCGACATGAAAGGCGACCTGCAGAATGACGACGGCACTGTAGTAGAGGTTTCCTCTATTCACCGTGCAGCACCTTCGTATGAAGATCAGTCTAACTCGCAGGAACTGCTGGAAACCGGCATCAAGGTTATCGACCTGATGTGTCCGTTCGCTAAGGGCGGTAAAGTCGGTCTGTTCGGTGGTGCGGGTGTAGGTAAAACCGTAAACATGATGGAGCTGATCCGTAACATCGCGGCTGAGCACTCAGGTTATTCGGTATTTGCTGGTGTGGGTGAGCGTACTCGTGAGGGTAACGACTTCTACCACGAAATGACTGACTCCAACGTTATCGAT---------------------AAAGTAGCGCTGGTGTATGGCCAGATGAACGAGCCGCCGGGTAACCGTCTGCGCGTAGCACTGACCGGTTTGACCATGGCGGAAAAATTCCGTGATGAAGGCCGTGACGTTCTGCTGTTCATCGACAACATCTATCGTTACACCCTGGCCGGTACTGAAGTATCAGCACTGCTGGGTCGTATGCCATCTGCGGTAGGCTATCAGCCAACGCTGGCAGAAGAGATGGGTGTGCGCATTCCAACACTGGAAGAGTGCGATGTCTGCCACGGTAGCGGCGCGAAAGCGGGGACCAAACCGCAGACCTGTCATACCTGTCATGGCGCAGGCCAGGTGCAGATGCGTCAGGGCTTCTTCACTGTGCAGCAGGCGTGTCCGACCTGTCACGGTCGCGGTTCAGTGATCAAAGATCCGTGCAATGCTTGTCATGGTCACGGTCGCGTTGAGCGCAGTAAAACCCTGTCGGTGAAAATTCCAGCAGGCGTGGATACCGGCGATCGCATTCGTCTGACCGGCGAAGGTGAAGCGGGCGAACAGGGCGCACCAGCGGGCGATCTGTACGTTCAGGTTTCGGTGAAAAAGCACCCGATCTTTGAGCGTGAAGATAACAACCTATATTGCGAAGTGCCGATTAACTTTGCGATGGCAGCATTGGGTGGCGAGATTGAAGTGCCGACGCTTGATGGGCGTGTGAACCTGAAAGTGCCTTCTGAAACGCAAACTGGCAAGCTGTTCCGCATGCGCGGTAAAGGCGTGAAATCGGTGCGTGGTGGTGCGGTAGGCGATTTGCTGTGTCGTGTGGTGGTGGAAACGCCAGTTAGCCTCAATGACAAACAGAAAGCGTTACTGCGTGAACTGGAAGAGAGTTTTGGCGGCCCGAGCGGTGAGAAAAACGTCGTAAACGCCCTGTCACAGAAGCTGGAGCTGACCATTCGCCGTGAAGGCAAAGTGCATCAGCAGGTTTATCAGCACGGCGTGCCGCAGGCACCGCTGGCGGTGTCCGGTGATACCGATGCAACCGGTACTCGCGTGCGTTTCTGGCCGAGCTACGAAACCTTCACCAATGTGATTGAGTTTGAGTACGAAATCCTGGCGAAACGCCTGCGTGAACTGTCGTTCCTGAACTCTGGCGTTTCGATTCGTCTGGAAGACAAACGCGAC---GGCAAGAACGATCACTTCCACTACGAAGGCGGCATCAAGGCGTTCGTTGAGTATCTCAACAAGAACAAAACCCCGATTCACCCAACGGTGTTCTACTTCTCGACGGAG---AAAGATGGCATTGGCGTGGAAGTGGCGCTGCAGTGGAACGATGGTTTCCAGGAAAACATCTACTGCTTCACCAACAACATTCCACAGCGCGACGGCGGTACGCACCTGGCGGGCTTCCGTGCGGCAATGACGCGTACGCTGAATGCCTACATGGATAAAGAAGGCTACAGCAAAAAAGCCAAA------GTCAGTGCGACCGGTGACGATGCGCGTGAAGGCCTGATTGCAGTGGTTTCCGTGAAAGTGCCGGATCCGAAATTCTCTTCTCAGACCAAAGATAAGCTGGTCTCTTCTGAAGTGAAATCGGCGGTTGAGCAGCAGATGAACGAACTGCTGGCGGAATACCTGCTGGAAAATCCGTCTGACGCGAAAATC
""")
    apath = "/Users/brett/tmp/test.phy"
    output_path = "/Users/brett/tmp/output"
    alignment.write(apath)
    # a = Analysis(apath, output_path, False)
    a = Analysis(apath, output_path, True, threads=-1)
    # a = Analysis(apath, output_path, True)

    pa = Partition('a', (1, 1000, 3))
    pb = Partition('b', (2, 1000, 3))
    pc = Partition('c', (3, 1000, 3))
    # s1 = Subset(pa, pb)
    # s2 = Subset(pc)
    # sch = Scheme('a', s1, s2)

    models = "JC+I JC K80 TrNef K81".split()
    # models = phyml_models.get_all_models()
    # a.analyse_all_schemes(models)
    a.analyse_all_possible(models)


