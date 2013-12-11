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

import logging
log = logging.getLogger("reporter")

import os
import itertools
import operator

scheme_header_template = "%-18s: %s\n"
scheme_subset_template = "%-6s | %-10s | %-10s | %-100s\n"
subset_template = "%-15s | %-15s | %-15s | %-15s | %-15s\n"


class TextReporter(object):
    def __init__(self, config):
        self.cfg = config
        self.cfg.reporter = self

    def write_subset_summary(self, sub):
        pth = os.path.join(self.cfg.subsets_path, sub.name + '.txt')
        # Sort everything
        model_results = [(r.bic, r) for r in sub.results.values()]
        model_results.sort()
        output = open(pth, 'w')
        # TODO change back to full name...
        # output.write("Model selection results for subset: %s\n" % sub.full_name)
        output.write("Model selection results for subset: %s\n" % sub.name)
        output.write("Subset alignment stored here: %s\n" % sub.alignment_path)
        output.write("This subset contains the following data_blocks: %s\n" % sub)
        output.write("Models are organised according to their BIC scores\n\n")
        output.write(subset_template % ("Model", "lNL", "AIC", "AICc", "BIC"))
        for bic, r in model_results:
            output.write(subset_template % (r.model, r.lnl, r.aic, r.aicc, r.bic))

    def write_scheme_summary(self, sch, result):
        pth = os.path.join(self.cfg.schemes_path, sch.name + '.txt')
        output = open(pth, 'w')
        self.output_scheme(sch, result, output)

    def output_scheme(self, sch, result, output):
        self.write_scheme_header(sch, result, output)
        sorted_subsets = [sub for sub in sch]
        sorted_subsets.sort(key=lambda sub: min(sub.columns), reverse=False)
        self.write_subsets(sch, result, output, sorted_subsets)
        self.write_nexus_summary(output, sorted_subsets)
        self.write_raxml(sch, result, output, sorted_subsets)

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
            partition_sites = sub.site_description

            output.write("\tcharset Subset%s = %s;\n" % (subset_number, partition_sites))
            charpartition.append("Group%s:Subset%s" % (subset_number, subset_number))
            subset_number += 1
        output.write('\tcharpartition PartitionFinder = %s;\n' % ', '.join(charpartition))
        output.write('end;\n')

    def write_subsets(self, sch, result, output, sorted_subsets):
        output.write(scheme_subset_template % (
            "Subset", "Best Model", "# sites", "Partition names"))
        number = 1
        # a way to print out the scheme in PF format
        pf_scheme_description = []
        

        for sub in sorted_subsets:
            partition_names = sub.long_name
            pf_scheme_description.append("(%s)" % partition_names)
            partition_sites = sub.site_description
            output.write(scheme_subset_template % (
                number, 
                sub.best_model, 
                len(sub.columns), 
                partition_names,
                ))
            number += 1

        pf_scheme_description = " ".join(pf_scheme_description)
        output.write("\n\nScheme Description in PartitionFinder format\n")
        output.write("Scheme_%s = %s;" % (sch.name, pf_scheme_description))

    def write_raxml(self, sch, result, output, sorted_subsets):
        """Print out partition definitions in RaxML-like format, might be
        useful to some people
        """
        from raxml_models import get_raxml_protein_modelstring
        output.write("\n\nRaxML-style partition definitions\n")

        subset_number = 1
        for sub in sorted_subsets:
            partition_sites = sub.site_description

            if self.cfg.datatype == "DNA":
                model = 'DNA'
            elif self.cfg.datatype == "protein":
                model = get_raxml_protein_modelstring(each_s.best_model)
            else:
                raise RuntimeError

            output.write("%s, Subset%s = %s" % (model, subset_number, partition_sites))
            output.write("\n")
            subset_number += 1

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
        if self.cfg.search in ["rcluster", "hcluster"]:
            pretty_weights = "rate = %s, base = %s, model = %s, alpha = %s" %(
                               str(self.cfg.cluster_weights["rate"]),
                               str(self.cfg.cluster_weights["freqs"]),
                               str(self.cfg.cluster_weights["model"]),
                               str(self.cfg.cluster_weights["alpha"]))
            output.write(scheme_header_template % ("weights", pretty_weights))
        if self.cfg.search == "rcluster":
            output.write(scheme_header_template % ("rcluster-percent",
                                                   self.cfg.cluster_percent))
            if self.cfg.cluster_max != None:
                output.write(scheme_header_template % ("rcluster-max",
                                                       self.cfg.cluster_max))
            else:
                output.write(scheme_header_template % ("rcluster-max",
                                                       "all rcluster-percent schemes"))

        output.write('\n\nBest partitioning scheme\n\n')
        self.output_scheme(result.best_scheme, result.best_result, output)
        log.info("Information on best scheme is here: %s", pth)
