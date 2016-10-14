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

import os
import shutil
from database import Database

from alignment import Alignment, SubsetAlignment
import threadpool
import scheme
import subset_ops
import results
import threading
import collections
from config import the_config
from util import PartitionFinderError, ExternalProgramError
import util
import raxml
from shutil import copyfile

class AnalysisError(PartitionFinderError):
    pass


class Analysis(object):
    """Performs the analysis and collects the results"""
    def __init__(self, cfg, force_restart, threads):
        the_config.validate()

        # TODO: Remove -- put this all into "options"
        if threads == -1:
            threads = threadpool.get_cpu_count()

        self.threads = threads

        # TODO: Move these to the config validate and prepare
        log.info("Beginning Analysis")
        self.process_restart(force_restart)

        # Make some folders for the analysis
        the_config.make_output_folders()
        the_config.database = Database(the_config)

        # Check for old analyses to see if we can use the old data
        the_config.check_for_old_config()

        # TODO: This is going to be in "Prepare"
        self.make_alignment(cfg.alignment_path)
        self.make_tree(cfg.user_tree_topology_path)

        # We need this to block the threads for critical stuff
        self.lock = threading.Condition(threading.Lock())

        # Store the result in here
        self.results = results.AnalysisResults(the_config.model_selection)

    def process_restart(self, force_restart):
        if force_restart:
            # Remove everything
            if os.path.exists(the_config.output_path):
                log.warning("Deleting all previous workings in '%s'" %
                            the_config.output_path)
                shutil.rmtree(the_config.output_path)
        else:
            # Remove the schemes folder, and clean out the phylofiles folder
            if os.path.exists(the_config.schemes_path):
                log.debug("Removing files in '%s'" % the_config.schemes_path)
                shutil.rmtree(the_config.schemes_path)
            if os.path.exists(the_config.phylofiles_path):
                log.debug("Removing files in '%s'" % the_config.phylofiles_path)
                shutil.rmtree(the_config.phylofiles_path)


    def analyse(self):
        try:
            self.do_analysis()
        finally:
            # TODO: Not really the right place for it?
            the_config.database.close()
        return self.results



    def make_alignment(self, source_alignment_path):
        # Make the alignment
        self.alignment = Alignment()
        self.alignment.read(source_alignment_path)

        # TODO REMOVE -- this should be part of the checking procedure
        # We start by copying the alignment
        self.alignment_path = os.path.join(the_config.start_tree_path, 'source.phy')
        if os.path.exists(self.alignment_path):
            # Make sure it is the same
            old_align = Alignment()
            old_align.read(self.alignment_path)
            if not old_align.same_as(self.alignment):
                log.error("""Alignment file has changed since previous run. You
                     need to use the force-restart option.""")
                raise AnalysisError

            compare = lambda x, y: collections.Counter(x) == collections.Counter(y)

            if not compare(old_align.species, self.alignment.species):
                log.error("""Species names in alignment have changed since previous run. You
                     need to use the force-restart option.""")
                raise AnalysisError


        else:
            self.alignment.write(self.alignment_path)

    def need_new_tree(self, tree_path):
        if os.path.exists(tree_path):
            if ';' in open(tree_path).read():
                log.info("Starting tree file found.")
                redo_tree = False
            else:
                log.info("""Starting tree file found but it is incomplete.
                             Re-estimating""")
                redo_tree = True
        else:
            log.info("Starting tree will be estimated from the data.")
            redo_tree = True

        return redo_tree

    def make_tree(self, user_path):
        # Begin by making a filtered alignment, containing ONLY those columns
        # that are defined in the subsets
        subset_with_everything = subset_ops.merge_subsets(the_config.user_subsets)
        self.filtered_alignment = SubsetAlignment(
            self.alignment, subset_with_everything)
        self.filtered_alignment_path = os.path.join(
            the_config.start_tree_path,  'filtered_source.phy')
        self.filtered_alignment.write(self.filtered_alignment_path)

        # Check the full subset against the alignment
        subset_ops.check_against_alignment(subset_with_everything, self.alignment, the_config)

        # We start by copying the alignment
        self.alignment_path = os.path.join(
            the_config.start_tree_path, 'source.phy')

        # Now check for the tree
        tree_path = the_config.processor.make_tree_path(
            self.filtered_alignment_path)

        if self.need_new_tree(tree_path):
            log.debug("Estimating new starting tree, no old tree found")

            # If we have a user tree, then use that, otherwise, create a topology
            util.clean_out_folder(the_config.start_tree_path,
                                  keep=["filtered_source.phy", "source.phy"])

            if user_path is not None and user_path != "":
                # Copy it into the start tree folder
                log.info("Using user supplied topology at %s" % user_path)
                topology_path = os.path.join(the_config.start_tree_path, 'user_topology.phy')
                util.dupfile(user_path, topology_path)
                need_bl = True
            elif the_config.no_ml_tree == True:
                log.debug(
                    "didn't find tree at %s, making a new one" % tree_path)
                topology_path = the_config.processor.make_topology(
                    self.filtered_alignment_path, the_config.datatype, the_config.cmdline_extras)
                need_bl = True
            elif the_config.no_ml_tree == False:
                log.debug(
                    "didn't find tree at %s, making an ML tree with RAxML" % tree_path)

                tree_scheme = scheme.create_scheme(
                    the_config, "tree_scheme", range(len(the_config.user_subsets)))

                topology_path = raxml.make_ml_topology(
                    self.filtered_alignment_path, the_config.datatype, the_config.cmdline_extras, tree_scheme, self.threads)
                
                # here we copy the ML tree topology so it can be used with PhyML too
                # TODO: this is a hack, and it would be better to decide on a universal
                # name for the different types of tree we might have.
                phyml_tree = os.path.join(os.path.dirname(topology_path), "filtered_source.phy_phyml_tree.txt")
                copyfile(topology_path, phyml_tree)

                need_bl = False

            if need_bl == True:
                # Now estimate branch lengths
                tree_path = the_config.processor.make_branch_lengths(
                    self.filtered_alignment_path,
                    topology_path,
                    the_config.datatype,
                    the_config.cmdline_extras)

        self.tree_path = tree_path
        log.debug("Starting tree with branch lengths is here: %s" %
                 self.tree_path)

    def run_task(self, model_name, sub):
        # This bit should run in parallel (forking the processor)
        try:
            the_config.processor.analyse(
                model_name,
                sub.alignment_path,
                self.tree_path,
                the_config.branchlengths,
                the_config.cmdline_extras
            )
            fabricate = False
        except ExternalProgramError:
            if not the_config.suppress_errors:
                # In the Kmeans algorithm we suppress errors and "fabricate"
                # subsets (we assume the error is because the subset is too
                # small for analysis)
                raise

            # If it is kmeans we assume that the error is because the subset
            # is too small or unanalysable, so we fabricate it
            log.debug("New subset could not be analysed. It will be merged "
                        "at the end of the analysis")
            fabricate = True

        # Not entirely sure that WE NEED to block here, but it is safer to do
        # It shouldn't hold things up toooo long...
        self.lock.acquire()
        try:
            if fabricate:
                sub.fabricate_model_result(the_config, model_name)
            else:
                sub.parse_model_result(the_config, model_name)

            # Try finalising, then the result will get written out earlier...
            sub.finalise(the_config)
        finally:
            self.lock.release()

    def add_tasks_for_sub(self, tasks, sub):
        for m in sub.models_to_process:
            tasks.append((self.run_task, (m, sub)))

    def run_concurrent(self, tasks):
        for func, args in tasks:
            func(*args)

    def run_threaded(self, tasks):
        if not tasks:
            return
        pool = threadpool.Pool(tasks, self.threads)
        pool.join()

    def analyse_list_of_subsets(self, all_subsets, ):
        # get a whole list of subsets analysed in parallel

        # analyse bigger subsets first, for efficiency
        all_subsets.sort(key = lambda x: 1.0/float(len(x.columns)))

        # chunk the list into blocks of ~1000 tasks 
        # in empirical testing, this speeds things up lot
        # though we are not entirely sure why...
        n = 1000        
        n = int(n / len(the_config.models))
        if(n<1): n=1 # seems unlikely...

        log.debug("chunk size (in number of subsets) = %d", n)

        subset_chunks = [all_subsets[i:i + n] for i in xrange(0, len(all_subsets), n)]
        
        for subsets in subset_chunks:
            # prepare the list of tasks
            tasks = []
            for sub in subsets:
                if sub.is_done:
                    pass
                elif sub.is_prepared:
                    self.add_tasks_for_sub(tasks, sub)
                else:
                    sub.prepare(the_config, self.alignment)
                    self.add_tasks_for_sub(tasks, sub)
            if tasks:
                # Now do the analysis
                if self.threads == 1:
                    self.run_concurrent(tasks)
                else:
                    self.run_threaded(tasks)

        # Now see if we're done
        for sub in all_subsets:
            # ALL subsets should already be finalised in the task. We just
            # check again here
            if not sub.finalise(the_config):
                log.error("Failed to run models %s; not sure why" %
                          ", " "".join(list(sub.models_not_done)))
                raise AnalysisError

    def analyse_scheme(self, sch):
        # Progress
        the_config.progress.next_scheme()

        # analyse the subsets in the scheme that aren't done
        # NB for most schemes we will have all subsets done, so this saves time
        not_done = []
        for sub in sch:
            if sub.is_done == False:
                not_done.append(sub)
        if not_done:
            self.analyse_list_of_subsets(not_done)

        # AIC needs the number of sequences
        number_of_seq = len(self.alignment.species)
        result = scheme.SchemeResult(sch, number_of_seq, the_config.branchlengths, the_config.model_selection)
        self.results.add_scheme_result(sch, result)

        return result
