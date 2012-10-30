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

"""This is where everything comes together, and we do the analysis"""

import logging
log = logging.getLogger("analysis")

import os, shutil

from alignment import Alignment, SubsetAlignment
import phyml, phyml_models
import raxml, raxml_models
import threadpool
import scheme
import subset
import util
import results

from util import PartitionFinderError
class AnalysisError(PartitionFinderError):
    pass

class Analysis(object):
    """Performs the analysis and collects the results"""
    def __init__(self, cfg, rpt, 
                 force_restart=False, 
                 save_phylofiles=False,
                 phylogeny_program = 'phyml',
                 threads=-1):
        cfg.validate()

        self.cfg = cfg
        self.rpt = rpt
        self.threads = threads
        self.save_phylofiles = save_phylofiles
        self.results = results.AnalysisResults()
    

        #TODO this is very ugly, it would liek to be prettier
        if phylogeny_program=='phyml' or phylogeny_program==-1:
            self.processor = phyml
            self.models = phyml_models
            self.get_model_difficulty = phyml_models.get_model_difficulty
        elif phylogeny_program=='raxml':
            self.processor = raxml
            self.models = raxml_models
            self.get_model_difficulty = raxml_models.get_model_difficulty
        else:
            log.error("Unrecognised option %s for phylogeny program, only PhyML and RAxML are "
                      "currently supported" % phylogeny_program)
            raise AnalysisError
            
        log.info("Beginning Analysis")
        if force_restart:
            # Remove everything
            if os.path.exists(self.cfg.output_path):
                log.warning("Deleting all previous workings in '%s'", 
                            self.cfg.output_path)
                shutil.rmtree(self.cfg.output_path)
        else:
            # Just remove the schemes folder
            if os.path.exists(self.cfg.schemes_path):
                log.info("Removing Schemes in '%s' (they will be "
                         "recalculated from existing subset data)",
                         self.cfg.schemes_path)
                shutil.rmtree(self.cfg.schemes_path)

        #check for old analyses to see if we can use the old data
        self.cfg.check_for_old_config()

        # Make some folders for the analysis
        self.cfg.make_output_folders()
        self.make_alignment(cfg.alignment_path)

        self.make_tree(cfg.user_tree_topology_path)
        self.subsets_analysed_set = set() #a counter for user info
        self.subsets_analysed = 0 #a counter for user info
        self.total_subset_num = None
        self.schemes_analysed = 0 #a counter for user info
        self.total_scheme_num = None

    def analyse(self):
        self.do_analysis()
        self.results.finalise()
        self.report()
        return self.results

    def report(self):
        best = [
            ("Best scheme according to AIC", self.results.best_aic),
            ("Best scheme according to AICc", self.results.best_aicc),
            ("Best scheme according to BIC", self.results.best_bic),
        ]
        self.rpt.write_best_schemes(best)
        self.rpt.write_all_schemes(self.results)

    def make_alignment(self, source_alignment_path):
        # Make the alignment 
        self.alignment = Alignment()
        self.alignment.read(source_alignment_path)

        # We start by copying the alignment
        self.alignment_path = os.path.join(self.cfg.start_tree_path, 'source.phy')
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

    def make_tree(self, user_path):
        # Begin by making a filtered alignment, containing ONLY those columns
        # that are defined in the subsets
        subset_with_everything = subset.Subset(*list(self.cfg.partitions))
        self.filtered_alignment = SubsetAlignment(self.alignment, 
                                                  subset_with_everything)
        self.filtered_alignment_path = os.path.join(self.cfg.start_tree_path,
                                                    'filtered_source.phy')
        self.filtered_alignment.write(self.filtered_alignment_path)

        # Now we've written this alignment, we need to lock everything in
        # place, no more adding partitions, or changing them from now on.
        self.cfg.partitions.check_against_alignment(self.alignment)
        self.cfg.partitions.finalise()

        # We start by copying the alignment
        self.alignment_path = os.path.join(self.cfg.start_tree_path, 'source.phy')

        # Now check for the tree
        tree_path = self.processor.make_tree_path(self.filtered_alignment_path)
        if not os.path.exists(tree_path):
            # If we have a user tree, then use that, otherwise, create a topology
            if user_path != None and user_path != "":
                # Copy it into the start tree folder
                log.info("Using user supplied topology at %s", user_path)
                topology_path = os.path.join(self.cfg.start_tree_path,
                                             'user_topology.phy')
                self.processor.dupfile(user_path, topology_path)
            else:
                topology_path = self.processor.make_topology(self.filtered_alignment_path, self.cfg.datatype)

            # Now estimate branch lengths
            tree_path = self.processor.make_branch_lengths(self.filtered_alignment_path, topology_path, self.cfg.datatype)
                
        self.tree_path = tree_path
        log.info("Starting tree with branch lengths is here: %s", self.tree_path) 


    def analyse_subset(self, sub, models):
        """Analyse the subset using the models given
        This is the core place where everything comes together
        The results are placed into subset.result
        """

        log.debug("About to analyse %s using models %s", sub, ", ".join(list(models)))

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

        subset_cache_path = os.path.join(self.cfg.subsets_path, sub.name + '.bin')
        # We might have already saved a bunch of results, try there first
        if not sub.results:
            log.debug("Reading in cached data from the subsets file")
            sub.read_cache(subset_cache_path)

        # First, see if we've already got the results loaded. Then we can
        # shortcut all the other checks
        models_done = set(sub.results.keys())
        log.debug("These models have already been done: %s", models_done)
        models_required = set(models)
        models_to_do = models_required - models_done
        log.debug("Which leaves these models still to analyse: %s", models_to_do)

        

        
        # Empty set means we're done
        if not models_to_do:
            log.debug("All models already done, so using just the cached results for subset %s", sub)
            #if models_done!=set(models): #redo model selection if we have different models
            sub.model_selection(self.cfg.model_selection, self.cfg.models)        
            return


        # Make an Alignment from the source, using this subset
        sub_alignment = SubsetAlignment(self.alignment, sub)
        sub_path = os.path.join(self.cfg.phyml_path, sub.name + '.phy')
        # Add it into the sub, so we keep it around
        sub.alignment_path = sub_path

        # Maybe it is there already?
        if os.path.exists(sub_path):
            log.debug("Found existing alignment file %s", sub_path)
            old_align = Alignment()
            old_align.read(sub_path)

            # It had better be the same!
            if not old_align.same_as(sub_alignment):
                log.error("It looks like you have changed one or more of the"
                        "data_blocks in the configuration file, "
                        "so the new subset alignments"
                        " don't match the ones stored for this analysis."
                        "You'll need to run the program with --force-restart")
                raise AnalysisError
        else:
            # We need to write it
            sub_alignment.write(sub_path)

        # Try and read in some previous analyses
        log.debug("Checking for old results in the phyml folder")
        self.parse_results(sub, models_to_do)
        if not models_to_do:
            #if models_done!=set(models): #redo model selection if we have different models
            sub.model_selection(self.cfg.model_selection, self.cfg.models)        
            return

        # What is left, we actually have to analyse...
        tasks = []

        #for efficiency, we rank the models by their difficulty - most difficult first
        difficulty = []        
        for m in models_to_do:
            difficulty.append(self.get_model_difficulty(m))
        
        #hat tip to http://scienceoss.com/sort-one-list-by-another-list/
        difficulty_and_m = zip(difficulty, models_to_do)
        difficulty_and_m.sort(reverse=True)
        sorted_difficulty, sorted_models_to_do = zip(*difficulty_and_m)
            
        log.debug("About to analyse these models, in this order: %s", sorted_models_to_do)
        for m in sorted_models_to_do:
            #a_path, out_path = phyml.make_analysis_path(self.cfg.phyml_path, sub.name, m)
            tasks.append((self.processor.analyse, 
                          (m, sub_path, self.tree_path, self.cfg.branchlengths)))

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
        # but ONLY on the models specified in the cfg file.
        sub.model_selection(self.cfg.model_selection, self.cfg.models)        
        
        # If we made it to here, we should write out the new summary
        self.rpt.write_subset_summary(sub)
        # We also need to update this
        sub.write_cache(subset_cache_path)

    def parse_results(self, sub, models_to_do):
        """Read in the results and parse them"""
        models_done = []
        for m in list(models_to_do):
            # sub.alignment_path
            stats_path, tree_path = self.processor.make_output_path(sub.alignment_path, m)
            if os.path.exists(stats_path):
                sub_output = open(stats_path, 'rb').read()
                # Annotate with the parameters of the model
                try:
                    result = self.processor.parse(sub_output)
                    sub.add_model_result(m, result, self.models)
                    # Remove the current model from remaining ones
                    models_to_do.remove(m)
                    
                    # Just used for below
                    models_done.append(m)
                    if self.save_phylofiles:
                        pass
                    else:
                        #we have the ifs because raxml doesn't output a treefile in some cases
                        if os.path.isfile(stats_path):
                            os.remove(stats_path)
                        if os.path.isfile(tree_path):
                            os.remove(tree_path)

                except self.processor.PhylogenyProgramError:
                    log.warning("Failed loading parse output from %s."
                              "Output maybe corrupted. I'll run it again.",
                              stats_path)

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
        result = scheme.SchemeResult(sch, number_of_seq, self.cfg.branchlengths)
        self.results.add_scheme_result(result)

        # TODO: should put all paths into config. Then reporter should decide
        # whether to create stuff
        fname = os.path.join(self.cfg.schemes_path, sch.name+'.txt')
        self.rpt.write_scheme_summary(result, open(fname, 'w'))

        return result

