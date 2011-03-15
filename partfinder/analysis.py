"""This is where everything comes together, and we do the analysis"""

import logging
log = logging.getLogger("analysis")

import os

from alignment import Alignment, SubsetAlignment
import scheme, subset, partition, phyml
from config import settings

class AnalysisError(Exception):
    pass

def _make_folder(pth):
    if os.path.exists(pth):
        if not os.path.isdir(pth):
            log.error("Cannot create folder '%s'", pth)
            raise AnalysisError
    else:
        os.mkdir(pth)

# TODO HM. maybe don't need to initialise -- can do it all here?
class Analysis(object):
    """Performs the analysis and collects the results"""
    def __init__(self, alignment_path, output_path, force_restart):

        # Ignore the force restart for now
        self.output_path = output_path
        _make_folder(self.output_path)

        # Now make the alignment and tree
        self.alignment = Alignment()
        self.alignment.read(alignment_path)
        # We start by copying the alignment
        self.alignment_path = os.path.join(self.output_path, 'source.phy')
        if os.path.exists(self.alignment_path):
            # Make sure it is the same
            old_align = Alignment()
            old_align.read(self.alignment_path)
            if old_align != self.alignment:
                log.error("Alignment file has changed since previous run. "
                          "You need to use the force-restart option.")
                raise AnalysisError
        else:
            self.alignment.write(self.alignment_path)

        # Now check for the tree
        tree_path = phyml.tree_path(self.alignment_path)
        if not os.path.exists(tree_path):
            phyml.make_tree(self.alignment_path)

    def analyse_subset(sub):
        # Has it already been done?
        if sub.has_analysis:
            return

        # Make an Alignment from the source, using this subset
        sub_alignment = alignment.SubsetAlignment(settings.source, sub)

        # Let's get a filename
        alignment_pth = os.path.join(settings.output_path, sub.name + '.phy')

        # Maybe it is there already?
        if os.path.exists(alignment_pth):
            log.debug("Found existing alignment file %s", alignment_pth)
            old_align = SubsetAlignment(sub.name)
            old_align.read(alignment_pth)

            # It had better be the same!
            if not old_align.same_as(sub_alignment):
                log.error("It looks like you have changed something in the "
                        "configuration and I cannot trust the old analysis."
                        "You'll need to run the program with --force-restart")
                raise AnalysisError
        else:
            # We need to write it
            sub_alignment.write(alignment_pth)

        results = []
        for model in settings.models:
            pass
            # results.append(phyml.

        # Collect the subset results
        # Decide on the best and set it in the subset

    def analyse_scheme(sch):
        for sub in sch:
            analyse_subset(sub)

        # Now piece together the bits

    def analyse_all_schemes():
        """Process everything!"""
        # First make a tree
        #
        for sch in scheme.all_schemes:
            analyse_scheme(sch)

        # Now check the best
    

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    # tree_pth = phyml.make_tree(settings.program_path, settings.alignment)
    # settings.tree = tree_pth
    # settings.source = alignment.SourceAlignment(settings.alignment)

    # p1 = partition.Partition('one', (1, 10))
    # p2 = partition.Partition('two', (11, 20))
    # scheme.generate_all_schemes()
    # analyse_all_schemes()

