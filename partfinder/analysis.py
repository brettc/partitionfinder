#Copyright (C) 2012 Robert Lanfear and Brett Calcott
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
#program, the RAxML program, and the PyParsing library,
#all of which are protected by their own licenses and conditions, using
#PartitionFinder implies that you agree with those licences and conditions as well.

import logging
log = logging.getLogger("analysis")

import os
import shutil

from alignment import Alignment, SubsetAlignment
import threadpool
import scheme
import subset
import subset_ops
import results
import threading
from util import PartitionFinderError
import util
from util import PhylogenyProgramError

class AnalysisError(PartitionFinderError):
    pass


class Analysis(object):
    """Performs the analysis and collects the results"""
    def __init__(self, cfg, force_restart=False, threads=-1):
        cfg.validate()
        self.cfg = cfg
        self.threads = threads

        self.results = results.AnalysisResults(self.cfg.model_selection)
        log.info("Beginning Analysis")
        self.process_restart(force_restart)

        # Check for old analyses to see if we can use the old data
        self.cfg.check_for_old_config()
        # Make some folders for the analysis
        self.cfg.make_output_folders()
        self.make_alignment(cfg.alignment_path)
        self.make_tree(cfg.user_tree_topology_path)
        
        # We need this to block the threads for critical stuff
        self.lock = threading.Condition(threading.Lock())
    def process_restart(self, force_restart):
        if force_restart:
            # Remove everything
            if os.path.exists(self.cfg.output_path):
                log.warning("Deleting all previous workings in '%s'", self.cfg.output_path)
                shutil.rmtree(self.cfg.output_path)
        else:
            # Just remove the schemes folder
            if os.path.exists(self.cfg.schemes_path):
                log.info("Removing Schemes in '%s' (they will be recalculated from existing subset data)", self.cfg.schemes_path)
                shutil.rmtree(self.cfg.schemes_path)

    def analyse(self):
        self.do_analysis()
        return self.results

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
                log.error("Alignment file has changed since previous run. You need to use the force-restart option.")
                raise AnalysisError

        else:
            self.alignment.write(self.alignment_path)

    def need_new_tree(self, tree_path):
        if os.path.exists(tree_path):
            if ';' in open(tree_path).read():
                log.info("Starting tree file found.")
                redo_tree = False
            else:
                log.info("Starting tree file found but incomplete. Re-estimating")
                redo_tree = True
        else:
            log.info("Starting tree will be estimated from the data.")
            redo_tree = True

        return redo_tree

    def make_tree(self, user_path):
        # Begin by making a filtered alignment, containing ONLY those columns
        # that are defined in the subsets
        subset_with_everything = subset_ops.merge_subsets(self.cfg.user_subsets)
        self.filtered_alignment = SubsetAlignment(
            self.alignment, subset_with_everything)
        self.filtered_alignment_path = os.path.join(
            self.cfg.start_tree_path,  'filtered_source.phy')
        self.filtered_alignment.write(self.filtered_alignment_path)

        # TODO: This checking should still be done...
        # self.cfg.partitions.check_against_alignment(self.alignment)
        # self.cfg.partitions.finalise()

        # We start by copying the alignment
        self.alignment_path = os.path.join(
            self.cfg.start_tree_path, 'source.phy')

        # Now check for the tree
        tree_path = self.cfg.processor.make_tree_path(
            self.filtered_alignment_path)

        if self.need_new_tree(tree_path):
            log.debug("Estimating new starting tree, no old tree found")

            # If we have a user tree, then use that, otherwise, create a topology
            util.clean_out_folder(self.cfg.start_tree_path,
                                  keep=["filtered_source.phy", "source.phy"])

            if user_path is not None and user_path != "":
                # Copy it into the start tree folder
                log.info("Using user supplied topology at %s", user_path)
                topology_path = os.path.join(self.cfg.start_tree_path, 'user_topology.phy')
                self.cfg.processor.dupfile(user_path, topology_path)
            else:
                log.debug(
                    "didn't find tree at %s, making a new one" % tree_path)
                topology_path = self.cfg.processor.make_topology(
                    self.filtered_alignment_path, self.cfg.datatype, self.cfg.cmdline_extras)

            # Now estimate branch lengths
            tree_path = self.cfg.processor.make_branch_lengths(
                self.filtered_alignment_path,
                topology_path,
                self.cfg.datatype,
                self.cfg.cmdline_extras)

        self.tree_path = tree_path
        log.info("Starting tree with branch lengths is here: %s", self.tree_path)

    def run_task(self, model_name, sub):
        # This bit should run in parallel (forking the processor)
        try:
            self.cfg.processor.analyse(
                model_name,
                sub.alignment_path,
                self.tree_path,
                self.cfg.branchlengths,
                self.cfg.cmdline_extras
            )

        except PhylogenyProgramError as e:
            # TODO: probably should do something smart with these "errors" so
            # that we can pull them up and see what went wrong
            sub.analysis_error = e.stdout, e.stderr

        # Not entirely sure that WE NEED to block here, but it is safer to do
        # It shouldn't hold things up toooo long...
        self.lock.acquire()
        try:
            if sub.analysis_error == None:
                sub.parse_model_result(self.cfg, model_name)

            else:
                sub.fabricate_result(self.cfg, model_name)

            # Try finalising, then the result will get written out earlier...
            sub.finalise(self.cfg)
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

    def analyse_scheme(self, sch):
        # Progress
        self.cfg.progress.next_scheme()
    
        # Prepare by reading everything in first
        tasks = []
        for sub in sch:
            sub.prepare(self.cfg, self.alignment)
            self.add_tasks_for_sub(tasks, sub)
        # Now do the analysis
        if self.threads == 1:
            self.run_concurrent(tasks)
        else:
            self.run_threaded(tasks)

        # Now see if we're done
        for sub in sch:
            # ALL subsets should already be finalised in the task. We just
            # check again here
            if not sub.finalise(self.cfg):
                log.error("Failed to run models %s; not sure why", ", ".join(list(sub.models_to_do)))
                raise AnalysisError

        # AIC needs the number of sequences
        number_of_seq = len(self.alignment.species)
        result = scheme.SchemeResult(sch, number_of_seq, self.cfg.branchlengths, self.cfg.model_selection)
        self.results.add_scheme_result(sch, result)
        
        return result
