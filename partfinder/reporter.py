# Copyright (C) 2012 Robert Lanfear and Brett Calcott
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
import pandas
log = logtools.get_logger()
from config import the_config
from model_utils import *


import os

scheme_header_template = "%-18s: %s\n"
scheme_subset_template = "%-6s | %-10s | %-10s | %-32s | %-100s\n"
subset_template = "%-15s | %-15s | %-15s | %-15s  | %-15s | %-15s\n"

# We write different output for these searches
_odd_searches = ['kmeans', 'krmeans']
_scheme_data_csv = 'scheme_data.csv'

class TextReporter(object):
    def __init__(self, config):
        self.cfg = config
        self.cfg.reporter = self
        self.header_done = False


    def write_subset_summary(self, sub):
        pth = os.path.join(self.cfg.subsets_path, sub.subset_id + '.txt')
        # Sort everything

        cols = ['model_id', 'params', 'lnl', 'aicc', 'aic', 'bic']

        descr = self.cfg.data_layout.data_type.descr
        indices = dict([(t[0], i) for i, t in enumerate(descr) if t[0] in cols])

        sorted_results = [(row['aicc'], row) for row in sub.result_array]
        sorted_results.sort()


        output = open(pth, 'w')
        output.write("Model selection results for subset: %s\n" % sub.subset_id)
        if sub.alignment_path:
            output.write("Subset alignment stored here: %s\n" % sub.alignment_path)
        if the_config.search not in _odd_searches:
            output.write("This subset contains the following data_blocks: %s\n" % sub.name)
        output.write("Number of columns in subset: %d\n" % len(sub.columns))
        output.write("Models are organised according to their AICc scores\n\n")

        output.write(subset_template % ("Model", "Parameters", "lnL", "AICc", "AIC", "BIC"))
        for aicc, row in sorted_results:
            output.write(subset_template % (row[indices['model_id']], 
                                            row[indices['params']], 
                                            row[indices['lnl']], 
                                            row[indices['aicc']], 
                                            row[indices['aic']],
                                            row[indices['bic']]))  


    def write_scheme_summary(self, sch, result):
        pth = os.path.join(self.cfg.schemes_path, sch.name + '.txt')
        output = open(pth, 'w')
        self.output_scheme(sch, result, output)
        summary_pth = os.path.join(self.cfg.schemes_path, _scheme_data_csv)
        summary_output = open(summary_pth, "a")
        self.add_scheme_to_csv(sch, result, summary_output)

    def add_scheme_to_csv(self, sch, result, summary_output):
        if not self.header_done:
            summary_output.write('name,sites,lnL,parameters,subsets,aic,aicc,bic\n')
            self.header_done = True

        summary_output.write('%s,%d,%.2f,%d,%d,%.2f,%.2f,%.2f\n'
            %(sch.name,
                result.nsites,
                result.lnl,
                result.sum_k,
                result.nsubs,
                result.aic,
                result.aicc,
                result.bic
                )
            )

    def output_scheme(self, sch, result, output):
        self.write_scheme_header(sch, result, output)
        sorted_subsets = [sub for sub in sch]
        sorted_subsets.sort(key=lambda sub: min(sub.columns), reverse=False)
        self.write_subsets(sch, result, output, sorted_subsets)
        self.write_nexus_summary(output, sorted_subsets)
        self.write_raxml(sch, result, output, sorted_subsets)
        if self.cfg.datatype != "morphology":
            self.write_mrbayes(sch, result, output, sorted_subsets)

    def write_scheme_header(self, sch, result, output):
        output.write(scheme_header_template % ("Scheme Name", sch.name))
        output.write(scheme_header_template % ("Scheme lnL", result.lnl))
        if self.cfg.model_selection == "aic":
            output.write(scheme_header_template % ("Scheme AIC", result.aic))
        if self.cfg.model_selection == "aicc":
            output.write(scheme_header_template % ("Scheme AICc", result.aicc))
        if self.cfg.model_selection == "bic":
            output.write(scheme_header_template % ("Scheme BIC", result.bic))
        output.write(scheme_header_template % ("Number of params", result.sum_k))
        output.write(scheme_header_template % ("Number of sites", result.nsites))
        output.write(scheme_header_template % ("Number of subsets", result.nsubs))
        output.write("\n")

    def write_nexus_summary(self, output, sorted_subsets):
        output.write("\n\nNexus formatted character sets\n")
        output.write("begin sets;\n")

        subset_number = 1
        charpartition = []
        for sub in sorted_subsets:
            if self.cfg.search in _odd_searches:
                sites = [x + 1 for x in sub.columns]
                partition_sites = str(sites).strip('[]')
            else:
                partition_sites = sub.site_description_no_commas

            output.write("\tcharset Subset%s = %s;\n" % (subset_number, partition_sites))
            charpartition.append("Group%s:Subset%s" % (subset_number, subset_number))
            subset_number += 1
        output.write('\tcharpartition PartitionFinder = %s;\n' % ', '.join(charpartition))
        output.write('end;\n')

    def write_subsets(self, sch, result, output, sorted_subsets):
        
        output.write(scheme_subset_template % (
            "Subset", "Best Model", "# sites", "subset id", "Partition names"))
        number = 1
        # a way to print out the scheme in PF format
        pf_scheme_description = []
            
        if self.cfg.search in _odd_searches:
            for sub in sorted_subsets:
                num_sites = len(sub.columns)
                
                sites = [x + 1 for x in sub.columns]
                pf_scheme_description.append("(%s)" % str(sites).strip('[]'))
                output.write(scheme_subset_template % (
                    number, 
                    sub.best_model, 
                    num_sites,
                    sub.subset_id, 
                    'NA',
                    ))
                number += 1

        else:
            for sub in sorted_subsets:
                pf_scheme_description.append("(%s)" % sub.name)
                output.write(scheme_subset_template % (
                    number, 
                    sub.best_model, 
                    len(sub.columns), 
                    sub.subset_id, 
                    sub.name,
                    ))
                number += 1

        
        pf_scheme_description = " ".join(pf_scheme_description)
        output.write("\n\nScheme Description in PartitionFinder format\n")
        output.write("Scheme_%s = %s;" % (sch.name, pf_scheme_description))


    def write_raxml(self, sch, result, output, sorted_subsets):
        self.write_raxml_warning(output)
        write_raxml_partitions(sch, output, sorted_subsets)


    def write_raxml_warning(self, output):
        output.write("\n\nRaxML-style partition definitions\n")
        output.write("Warning: RAxML allows for only a single model of rate"
                     " heterogeneity in partitioned analyses. I.e. all "
                     "partitions must be assigned one of three types of model:" 
                     " No heterogeneity (e.g. GTR); +G (e.g. GTR+G); or "
                     "+I+G (e.g. GTR+I+G). If the best models for your dataset"
                     "contain different types of model for different subsets "
                     "you will need to decide on "
                     "the best rate heterogeneity model before you run "
                     "RAxML. If you prefer to do things more rigorously, "
                     "you can run separate PartitionFinder analyses for each type "
                     "of rate heterogenetity "
                     "Then choose the scheme with the lowest AIC/AICc/BIC score. "
                     "Note that these re-runs will be quick!\n\n" 
                    )

    def write_mrbayes(self, sch, result, output, sorted_subsets):
        """Print out partition definitions in MrBayes format, might be
        useful to some people
        """
        output.write("\n\nMrBayes block for partition definitions\n")

        output.write("Warning: MrBayes only allows a relatively small "
                     "collection of models. If any model in your analysis is not one that "
                     "is included in MrBayes (e.g. by setting nst = 1, 2, or "
                     "6 for DNA sequences; or is not in the available list of protein models for MrBayes)" 
                     "then this MrBayes block will just set that model "
                     "to nst = 6 for DNA, or 'wag' for Protein. Similarly, the only additional parameters "
                     "that this MrBayes block will include are +I and +G. Other "
                     " parameters, such as +F and +X, are ignored. " 
                     "If you want to use this MrBayes block for your analysis, "
                     "please make sure to check it carefully before you use it "
                     "we've done our best to make it accurate, but there may "
                     "be errors that remain!\n\n"
                    )

        output.write("begin mrbayes;\n\n")
        subset_number = 1
        charpartition = []
        for sub in sorted_subsets:
            if self.cfg.search in _odd_searches:
                sites = [x + 1 for x in sub.columns]
                partition_sites = str(sites).strip('[]').replace(',','')
            else:
                partition_sites = sub.site_description_no_commas

            output.write("\tcharset Subset%s = %s;\n" % (subset_number, partition_sites))
            charpartition.append("Subset%s" % (subset_number))
            subset_number += 1
        output.write('\n\tpartition PartitionFinder = %d:%s;\n' %(len(charpartition), ', '.join(charpartition)))
        output.write('\tset partition=PartitionFinder;\n\n')

        subset_number = 1
        for sub in sorted_subsets:

            if self.cfg.datatype == "DNA":
                model_text = get_mrbayes_modeltext_DNA(sub.best_model, subset_number)
            elif self.cfg.datatype == "protein":
                model_text = get_mrbayes_modeltext_protein(sub.best_model, subset_number)
            else:
                raise RuntimeError

            output.write(model_text)
            subset_number += 1

        if len(sorted_subsets)>1:
            output.write('\n\tprset applyto=(all) ratepr=variable;\n')
            output.write('\tunlink statefreq=(all) revmat=(all) shape=(all) pinvar=(all) tratio=(all);\n')

            if the_config.branchlengths == 'unlinked':
                output.write('\tunlink brlens=(all);\n')

        output.write('\nend;\n')



    def write_best_scheme(self, result):
        pth = os.path.join(self.cfg.output_path, 'best_scheme.txt')
        output = open(pth, 'wb')
        output.write('Settings used\n\n')
        output.write(scheme_header_template % ("alignment", self.cfg.alignment_path))
        output.write(scheme_header_template % ("branchlengths", self.cfg.branchlengths))
        output.write(scheme_header_template % ("models", ', '.join(self.cfg.models)))
        output.write(scheme_header_template % ("model_selection",
                                                self.cfg.model_selection))
        output.write(scheme_header_template % ("search", self.cfg.search))
        if self.cfg.search in ["rcluster", "hcluster", "rclusterf"]:
            pretty_weights = "rate = %s, base = %s, model = %s, alpha = %s" %(
                               str(self.cfg.cluster_weights["rate"]),
                               str(self.cfg.cluster_weights["freqs"]),
                               str(self.cfg.cluster_weights["model"]),
                               str(self.cfg.cluster_weights["alpha"]))
            output.write(scheme_header_template % ("weights", pretty_weights))

        if self.cfg.search.startswith("rcluster"):
            output.write(scheme_header_template % ("rcluster-percent",
                                                   self.cfg.cluster_percent))
            output.write(scheme_header_template % ("rcluster-max",
                                                       self.cfg.cluster_max))

        if self.cfg.search == "kmeans" or self.cfg.search == "krmeans":
            output.write(scheme_header_template % ("min-subset-size",
                                                   self.cfg.min_subset_size))
            output.write(scheme_header_template % ("kmeans based on",
                                                   self.cfg.kmeans))
            if self.cfg.all_states == True:
                output.write(scheme_header_template % ("--all_states setting used",
                                                       self.cfg.all_states))


        output.write('\n\nBest partitioning scheme\n\n')
        self.output_scheme(result.best_scheme, result.best_result, output)
        log.info("Information on best scheme is here: %s", pth)

        citation_text = write_citation_text(self)

        # now we write subset summaries for all the subsets in the best scheme
        for s in result.best_scheme:
            self.write_subset_summary(s)

        log.info("\n")
        log.info("\n")

        for c in citation_text:
            log.info("%s", c)
            output.write(c)


def write_raxml_partitions(sch, output, sorted_subsets, use_lg = False):
    """Print out partition definitions in RaxML-like format, might be
    useful to some people
    """

    subset_number = 1
    for sub in sorted_subsets:
        if the_config.search in _odd_searches:
            sites = [x + 1 for x in sub.columns]
            partition_sites = str(sites).strip('[]')
        else:
            partition_sites = sub.site_description

        if the_config.datatype == "DNA":
            model = 'DNA'
        elif the_config.datatype == "protein":
            if use_lg == False:
                model = get_raxml_protein_modelstring(sub.best_model)
            else:
                model = get_raxml_protein_modelstring("LG+G")
        elif the_config.datatype == "morphology":
            model = get_raxml_morphology_modelstring(sub.best_model)

        else:
            raise RuntimeError

        output.write("%s, Subset%s = %s" % (model, subset_number, partition_sites))
        output.write("\n")
        subset_number += 1


def write_citation_text(self):
    """Tell users which papers to cite"""

    citation_text = []

    ref_PF2 = ("Lanfear, R., Frandsen, P. B., Wright, A. M., Senfeld, T., Calcott, B. (2016) "
             "PartitionFinder 2: new methods for selecting partitioned models of evolution for" 
             "molecular and morphological phylogenetic analyses. "
             "Molecular biology and evolution. DOI: dx.doi.org/10.1093/molbev/msw260")

    ref_PF1 = ("Lanfear, R., Calcott, B., Ho, S. Y., & Guindon, S. (2012). "
             "PartitionFinder: combined selection of partitioning schemes "
             "and substitution models for phylogenetic analyses. "
             "Molecular biology and evolution, 29(6), 1695-1701.")

    ref_rcluster = ("Lanfear, R., Calcott, B., Kainer, D., Mayer, C., "
                    "& Stamatakis, A. (2014). Selecting optimal "
                    "partitioning schemes for phylogenomic datasets. "
                    "BMC evolutionary biology, 14(1), 82.")

    ref_rclusterf = ref_rcluster

    ref_kmeans = ("Frandsen, P. B., Calcott, B., Mayer, C., & Lanfear, R. "
                  "(2015). Automatic selection of partitioning schemes for "
                  "phylogenetic analyses using iterative k-means clustering "
                  "of site rates. BMC Evolutionary Biology, 15(1), 13.")

    ref_phyml = ("Guindon, S., Dufayard, J. F., Lefort, V., Anisimova, M., "
                 "Hordijk, W., & Gascuel, O. (2010). New algorithms and "
                 "methods to estimate maximum-likelihood phylogenies: "
                 "assessing the performance of PhyML 3.0. "
                 "Systematic biology, 59(3), 307-321.")

    ref_raxml = ("Stamatakis, A. (2014). RAxML version 8: a tool for "
                 "phylogenetic analysis and post-analysis of large phylogenies. "
                 "Bioinformatics, 30(9), 1312-1313.")

    ref_morph = ("Lewis, P. O. (2001). A likelihood approach to estimating "
                 "phylogeny from discrete morphological character data. "
                 "Systematic biology, 50(6), 913-925.")

    citation_text.append("\n\n\n*Citations for this analysis*\n")
    citation_text.append("-----------------------------")

    citation_text.append("\n")

    citation_text.append("If you use this analysis in your published "
        "work, please cite "
        "the following papers on which your analysis relied.\n")

    citation_text.append("\n")
    citation_text.append("For the version of PartitionFinder you used, "
                         "please cite:\n")

    citation_text.append("%s\n" % ref_PF2)

    citation_text.append("\n")
    citation_text.append("For the %s algorithm you used, please cite:\n" 
                         % (self.cfg.search))

    if self.cfg.search == "rcluster" or self.cfg.search == "hcluster":
        citation_text.append("%s\n" % ref_rcluster)

    elif self.cfg.search == "kmeans" or self.cfg.search == "krmeans":
        citation_text.append("%s\n" % ref_kmeans)

    elif self.cfg.search == "greedy":
        citation_text.append("%s\n" % ref_PF1)

    elif self.cfg.search == "rclusterf":
        citation_text.append("%s\n" % ref_rclusterf)


    citation_text.append("\n")
    if self.cfg.phylogeny_program == 'phyml':
        citation_text.append("Your analysis also used PhyML, so please cite:\n")
        citation_text.append("%s\n" % ref_phyml)

    elif self.cfg.phylogeny_program == 'raxml':
        citation_text.append("Your analysis also used RAxML, so please cite:\n")
        citation_text.append("%s\n" % ref_raxml)
    citation_text.append("\n")

    if self.cfg.datatype == 'morphology':
        citation_text.append("For the model of morphological evolution you used, please cite:\n")
        citation_text.append("%s\n" % ref_morph)

    return citation_text
