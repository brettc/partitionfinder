"""This is where everything comes together, and we do the analysis"""

import logging
log = logging.getLogger("analysis")

import os, shutil

from partition import all_partitions

from alignment import Alignment, SubsetAlignment
import phyml
import threadpool
import scheme

class AnalysisError(Exception):
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
    def __init__(self, 
                 alignment_path, 
                 output_path, 
                 force_restart,
                 threads=1):

        log.info("Beginning Analysis")
        self.threads = threads
        if force_restart:
            if os.path.exists(output_path):
                log.warning("Deleting all previous workings in '%s'", output_path)
                shutil.rmtree(output_path)

        # Make some folders for the analysis
        self.output_path = output_path
        make_dir(self.output_path)
        self.make_output_dir('subsets')
        self.make_output_dir('schemes')
        self.make_output_dir('phyml')

        # Now make the alignment and tree
        self.alignment = Alignment()
        self.alignment.read(alignment_path)

        # We start by copying the alignment
        self.alignment_path = os.path.join(self.output_path, 'source.phy')
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
        log.info("I'm working from a copy of the tree that is here: %s",
                 self.alignment_path)

        # Now check for the tree
        tree_path = phyml.make_tree_path(self.alignment_path)
        if not os.path.exists(tree_path):
            tree_path = phyml.make_tree(self.alignment_path)
        self.tree_path = tree_path
        log.info("The tree we're using is here: %s", self.tree_path) 

    def make_output_dir(self, name):
        new_path = os.path.join(self.output_path, name)
        make_dir(new_path)
        setattr(self, name+"_path", new_path)

    def analyse_subset(self, sub, models):
        """Analyse the subset using the models given
        This is the core place where everything comes together
        The results are placed into subset.result
        """

        sub_bin_path = os.path.join(self.subsets_path, sub.name + '.bin')
        # We might have already saved a bunch of results, try there first
        if not obj.results:
            sub.read_summary_binary(sub_bin_path)

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
                          (m, sub_path, a_path, self.tree_path)))

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
        for sub in sch:
            self.analyse_subset(sub, models)
 
        # AIC needs the number of sequences 
        number_of_seq = len(self.alignment.species)
        sch.assemble_results(number_of_seq)
        sch.write_summary(os.path.join(self.schemes_path, sch.name+'.txt'))

        log.info("Scheme %s analysed. AIC=%.2f", sch.name, sch.aic)
		

    def analyse_current_schemes(self, models):
        """Process everything!"""
        current_schemes = [s for s in scheme.all_schemes]
        for s in current_schemes:
            self.analyse_scheme(s, models)
        self.write_best_scheme(current_schemes)

    def analyse_all_possible(self, models):
        gen_schemes = scheme.generate_all_schemes()
        for s in gen_schemes:
            self.analyse_scheme(s, models)
        self.write_best_scheme(gen_schemes)

    def write_best_scheme(self, list_of_schemes):
        # Which is the best?
        sorted_schemes = [(s.aic, s) for s in list_of_schemes]
        sorted_schemes.sort()
        best = sorted_schemes[0][1]
        best.write_summary(os.path.join(self.output_path, 'best_scheme.txt'))
        
    def greedy_analysis(self, models):
        #do stuff
        start_scheme = [0, 1, 2]
		
    
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


