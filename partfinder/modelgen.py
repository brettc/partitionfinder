import logging
log = logging.getLogger("modelgen")

import config
import subprocess, shlex

from pyparsing import (
    Word, Dict, OneOrMore, alphas, nums, Suppress, Optional, Group, stringEnd,
    delimitedList, pythonStyleComment, ParseException, line, lineno, col,
    Keyword, ParserElement, ParseException)

def run_modelgen(alignment_path, tdir, num_gamma_cat=4):
    mg = config.settings.modelgen_path
    command = "java -jar %s %s %d" % (mg, alignment_path, num_gamma_cat)
    p = subprocess.Popen(shlex.split(command),
                         shell=False,
                         cwd=tdir,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    return stdout

if __name__ == '__main__':
    import tempfile, os, shutil
    logging.basicConfig(level=logging.DEBUG)
    align = r""">spp1
CTTGAGGTTCAGAATGGTAATGAA------GTGCTGG
>spp2
CTTGAGGTACAAAATGGTAATGAG------AGCCTGG
>spp3
CTTGAGGTACAGAATAACAGCGAG------AAGCTGG
>spp4
CTCGAGGTGAAAAATGGTGATGCT------CGTCTGG
"""

    config.initialise("~/tmp")
    tdir = tempfile.mkdtemp()
    print tdir
    outname = os.path.join(config.settings.output_path, 'temp.fasta')
    print outname
    f = open(outname, 'wb')
    f.write(align)
    f.close()
    
    outp = run_modelgen(outname, tdir)
    open('output.txt', 'w').write(outp)


