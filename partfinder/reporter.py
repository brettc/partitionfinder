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

scheme_header_template = "%-15s: %s\n"
scheme_subset_template = "%-6s | %-10s | %-30s | %-30s | %-40s\n"
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

        
    def garli_sub_text(self, sub, number):
        """In garli, the subset models are specified in one file, and the charsets in a 
        separate nexus file. The model definitions look like this:

        [model1]
        datatype = nucleotide
        ratematrix = ( 0 1 2 2 3 4 )
        statefrequencies = estimate
        ratehetmodel = gamma
        numratecats = 4
        invariantsites = none
      
        The point of this function is to take a subset, and extract information for
        these models, and return it
        
        """

        models = {
            "JC"        :   ("( 0 0 0 0 0 0 )", "equal"),
            "K80"       :   ("( 0 1 0 0 1 0 )", "equal"),
            "TrNef"     :   ("( 0 1 0 0 2 0 )", "equal"),
            "K81"       :   ("( 0 1 2 2 1 0 )", "equal"),
            "TVMef"     :   ("( 0 1 2 3 1 4 )", "equal"),
            "TIMef"     :   ("( 0 1 2 2 3 0 )", "equal"),
            "SYM"       :   ("( 0 1 2 3 4 5 )", "equal"),
            "F81"       :   ("( 0 0 0 0 0 0 )", "estimate"),
            "HKY"       :   ("( 0 1 0 0 1 0 )", "estimate"),
            "TrN"       :   ("( 0 1 0 0 2 0 )", "estimate"),  
            "K81uf"     :   ("( 0 1 2 2 1 0 )", "estimate"),
            "TVM"       :   ("( 0 1 2 3 1 4 )", "estimate"),
            "TIM"       :   ("( 0 1 2 2 3 0 )", "estimate"),
            "GTR"       :   ("( 0 1 2 3 4 5 )", "estimate"),
            "Dayhoff"   :   ("dayhoff", "dayhoff"),
            "JTT"       :   ("jones", "jones"),
            "WAG"       :   ("wag", "wag"),
            "mtREV"     :   ("mtrev", "mtrev"),
            "MtMam"     :   ("mtmam", "mtmam")
        }    
        
        #now we build up the description one by one
        header = "[model" + str(number) + "]"
        data = "datatype = " + self.cfg.datatype

        #now the model itself
        elements = sub.best_model.split("+")
        model_name = elements[0]
 
        # some yoga to make sure we catch models that aren't in the list above
        # typically, these are protein models like LG

        try:
            ratemat = "ratematrix = " + models[model_name][0]
            statefreq = "statefrequencies = " + models[model_name][1]
            #state frequencies, with a special catch for protein models
            if len(elements)>1 and "F" in elements[1:] and self.cfg.datatype=="protein":
                inv = "statefrequencies = empirical"
            warning = 0
        except KeyError: 
            ratemat = "ratematrix = fixed"
            statefreq = "statefrequencies = fixed"
            warning = 1
                    
        if len(elements)>1 and "G" in elements[1:]:
            ratemod = "ratehetmodel = gamma\nnumratecats = 4"
        else:  
            ratemod = "ratehetmodel = none\nnumratecats = 1"
        if len(elements)>1 and "I" in elements[1:]:
            inv = "invariantsites = estimate"
        else:
            inv = "invariantsites = none"
        
        final = "\n".join([header, data, ratemat, statefreq, ratemod, inv, "\n"])
        
        return(final, warning)
    
    def write_garli(self, output, sorted_subsets):
        output.write("\n\nGARLI model definitions\n")
        output.write("Please double check for accuracy.\n")
        output.write("These can be pasted into the garli.conf file.\n\n")
        
        warning = 0
        for i, sub in enumerate(sorted_subsets):
            model, warning = self.garli_sub_text(sub, i+1)
            output.write(model)
            warning += warning
            
        if warning>0:
            warning = ("N.B. There are rate matrices specified for your subsets that are"
                      " not specified in GARLI. You will need to specify these in the"
                      " correct order in your nexus file if you are going to run GARLI."
                      " See https://www.nescent.org/wg_garli/" "Specifying_a_custom_amino_acid_rate_matrix for more information.")
            output.write(warning)



    def write_scheme_summary(self, sch, result):
        pth = os.path.join(self.cfg.schemes_path, sch.name + '.txt')
        output = open(pth, 'w')
        self.output_scheme(sch, result, output)

    def output_scheme(self, sch, result, output):
        self.write_scheme_header(sch, result, output)
        sorted_subsets = [sub for sub in sch]
        sorted_subsets.sort(key=lambda sub: min(sub.columns), reverse=False)
        self.write_subsets(sch, result, output, sorted_subsets)
        self.write_raxml(sch, result, output, sorted_subsets)
        self.write_garli(output, sorted_subsets)

    def write_scheme_header(self, sch, result, output):
        output.write(scheme_header_template % ("Scheme Name", sch.name))
        output.write(scheme_header_template % ("Scheme lnL", result.lnl))
        output.write(scheme_header_template % ("Scheme AIC", result.aic))
        output.write(scheme_header_template % ("Scheme AICc", result.aicc))
        output.write(scheme_header_template % ("Scheme BIC", result.bic))
        output.write(scheme_header_template % ("Num params", result.sum_k))
        output.write(scheme_header_template % ("Num sites", result.nsites))
        output.write(scheme_header_template % ("Num subsets", result.nsubs))
        output.write("\n")

    def write_subsets(self, sch, result, output, sorted_subsets):
        output.write(scheme_subset_template % (
            "Subset", "Best Model", "Subset Partitions", "Subset Sites", "Alignment"))
        number = 1

        pf_scheme_description = []
            # a way to print out the scheme in PF format

        for sub in sorted_subsets:
            desc = {}
            names = []
            for part in sub:
                names.append(part.name)
                for subpart in part.description:  # loop through each sub-part of the partition
                    desc[subpart[0]] = subpart

            #pretty print the sites in the scheme
            desc_starts = desc.keys()
            desc_starts.sort()
            parts = []
            for key in desc_starts:
                part = desc[key]
                if part[2] == 1:
                    text = "%s-%s" % (part[0], part[1])
                else:
                    text = "%s-%s\\%s" % tuple(part)
                parts.append(text)
            parts = ', '.join(parts)

            names.sort()
            names = ', '.join(names)

            pf_scheme_description.append("(%s)" % names)

            output.write(scheme_subset_template % (
                number, sub.best_model, names, parts, sub.alignment_path))
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
        number = 1
        for sub in sorted_subsets:

            desc = {}
            names = []
            for part in sub:
                names.append(part.name)
                for subpart in part.description:  # loop through each sub-part of the partition
                    desc[subpart[0]] = subpart

            # Pretty print the sites in the scheme
            desc_starts = desc.keys()
            desc_starts.sort()
            parts = []
            for key in desc_starts:
                part = desc[key]
                if part[2] == 1:
                    text = "%s-%s" % (part[0], part[1])
                else:
                    text = "%s-%s\\%s" % tuple(part)
                parts.append(text)
            parts = ', '.join(parts)

            if self.cfg.datatype == "DNA":
                model = "DNA"
            elif self.cfg.datatype == "protein":
                model = get_raxml_protein_modelstring(sub.best_model)
            else:
                raise RuntimeError

            line = "%s, p%s = %s\n" % (model, number, parts)
            output.write(line)

            number += 1

    def write_best_scheme(self, txt, result):
        pth = os.path.join(self.cfg.output_path, 'best_scheme.txt')
        output = open(pth, 'wb')
        output.write(txt)
        output.write('\n\n')
        self.output_scheme(result.best_scheme, result.best_result, output)
        log.info("Information on best scheme is here: %s", pth)
