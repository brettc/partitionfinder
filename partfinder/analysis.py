"""This is where everything comes together, and we do the analysis"""

import logging
log = logging.getLogger("analysis")

import os

import alignment, scheme, subset, partition, phyml
from config import settings

class AnalysisError(Exception):
    pass

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
    # logging.basicConfig()
    logging.basicConfig(level=logging.DEBUG)
    import config
    config.initialise("~/tmp", True)
    settings.alignment = "test.phy"
    tree_pth = phyml.make_tree(settings.program_path, settings.alignment)
    settings.tree = tree_pth
    settings.source = alignment.SourceAlignment(settings.alignment)

    p1 = partition.Partition('one', (1, 10))
    p2 = partition.Partition('two', (11, 20))
    scheme.generate_all_schemes()
    analyse_all_schemes()

class Subset:
    # From Subset
    def analyse(self):
        # Check first to see if we've got the results, otherwise calculate and
        # cache them.
        if self.partitions in results_cache:
            log.debug("Returning cached result for %s", self)
            return results_cache[self.partitions]

        log.debug("Calculating result for %s", self)
        result = self._really_analyse()
        results_cache[self.partitions] = result
        return result

    def _really_analyse(self):
        fname = self.make_filename()
        sa = alignment.SubsetAlignment(
            fname, config.settings.source_alignment, self)
        return sa.analyse()

class X:
    def source_exists(self):
        return os.path.exists(self.source_path)
    def analysis_exists(self):
        return os.path.exists(self.analysis_path)

    def same_as_saved(self):
        spec, slen = self.read_source()
        if spec == self.species:
            log.debug("Are same")
            return True
        else:
            log.warning("different")
            return False

    def analyse(self):
        if self.source_exists():
            # We're already written a file of this name
            log.debug("%s already exists at '%s'", self, self.source_path)
            if self.same_as_saved():
                log.debug("%s Same", self)
                same_as = True
                XXXXXXXX
        else:
            # Otherwise write it out
            self.write_source()

        fresh_analysis = True
        if self.analysis_exists():
            log.debug("Reading in previous analysis of %s", self)
            output = file(self.analysis_path, 'r').read()
            fresh_analysis = False
        else:
            output = modelgen.run(self.source_path)
            log.debug("Saving analysis output of %s to %s", self,
                      self.analysis_path)
            open(self.analysis_path, 'w').write(output)

        log.debug("Parsing ModelGenerator output for %s", self)
        result = modelgen.parse(output)

        if fresh_analysis:
            # Show them how long it took.
            log.debug("New analysis of %s took %d seconds", self, result.processing_time)
        return result
    @property
    def source_path(self):
        return self.path + "." + alignment_format

    @property
    def analysis_path(self):
        return self.path + ".out"

    @property
    def path(self):
        # Defer to something that get's defined differently in subclasses
        # And cache it.
        if not hasattr(self, "_path"):
            self._path = self.get_path()
        return self._path

    def get_path(self):
        # Don't use this class -- use one of the subclasses below
        raise NotImplemented
