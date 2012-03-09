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

import logging
log = logging.getLogger("reporter")


scheme_header_template = "%-15s: %s\n"
scheme_subset_template = "%-6s | %-10s | %-30s | %-30s | %-40s\n"

class TextReporter(object):
    def __init__(self):
        pass

    # TODO. REMOVE 'write' from these. They could be output to a display
    def write_scheme_summary(self, sch, results, output):
        self.write_scheme_header(sch, results, output)
        # self.write_subsets(s)
        # self.write_raxml(s)

    def write_scheme_header(self, sch, results, output):
        output.write(scheme_header_template % ("Scheme Name", sch.name))
        output.write(scheme_header_template % ("Scheme lnL", results.lnl))
        output.write(scheme_header_template % ("Scheme AIC", results.aic))
        output.write(scheme_header_template % ("Scheme AICc", results.aicc))
        output.write(scheme_header_template % ("Scheme BIC", results.bic))
        output.write(scheme_header_template % ("Num params", results.sum_k))
        output.write(scheme_header_template % ("Num sites", results.nsites))
        output.write(scheme_header_template % ("Num subsets", results.nsubs))
        output.write("\n")

    # TODO FIX THIS and add them on. Write extra lines into the calling thing?
    def write_subsets(self, s):
        s.write(scheme_subset_template % (
            "Subset", "Best Model", "Subset Partitions", "Subset Sites",  "Alignment"))
        number = 1
        sorted_subsets = [sub for sub in self]
        sorted_subsets.sort(key=lambda sub: min(sub.columns), reverse=False)
                
        for sub in sorted_subsets:
            desc = {}
            names= []
            for part in sub:
                desc[part.description[0][0]] = part.description[0] #dict keyed by first site in part
                names.append(part.name)

            #pretty print the sites in the scheme
            desc_starts = desc.keys()
            desc_starts.sort()            
            parts = []
            for key in desc_starts:
                part = desc[key]
                if part[2]==1:
                    text = "%s-%s" %(part[0], part[1])
                else:
                    text = "%s-%s\\%s" % tuple(part)
                parts.append(text)
            parts = ', '.join(parts)
            	
            names.sort()
            names = ', '.join(names)
			
            s.write(scheme_subset_template % (
                number, sub.best_model, names, parts, sub.alignment_path))
            number += 1

    def write_raxml(self, s):
        """Print out partition definitions in RaxML-like format, might be
        useful to some people
        """

        s.write("\n\nRaxML-style partition definitions\n")
        number = 1
        for sub in sorted_subsets:
            desc = {}
            names= []
            for part in sub:
                desc[part.description[0][0]] = part.description[0] #dict keyed by first site in part
                names.append(part.name)

            #pretty print the sites in the scheme
            desc_starts = desc.keys()
            desc_starts.sort()            
            parts = []
            for key in desc_starts:
                part = desc[key]
                if part[2]==1:
                    text = "%s-%s" %(part[0], part[1])
                else:
                    text = "%s-%s\\%s" % tuple(part)
                parts.append(text)
            parts = ', '.join(parts)
            line = "DNA, p%s = %s\n" %(number, parts)
            s.write(line)

            number += + 1
