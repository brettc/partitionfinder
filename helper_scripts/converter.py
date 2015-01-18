#!/usr/bin/env python
import dendropy
from dendropy import Tree, TaxonSet
from dendropy.utility.fileutils import find_files
import os
import sys

print sys.argv[1]
print sys.argv[2]
print sys.argv[3]
print sys.argv[4]
o_file = sys.argv[1]
flist = find_files(top=o_file, filename_filter=sys.argv[2])
bases = [os.path.splitext(filename)[0] for filename in flist]

datas = [dendropy.StandardCharacterMatrix.get_from_path(filename, sys.argv[3], preserve_underscores=True) for filename in flist]
[data.write_to_path(base, sys.argv[4], space_to_underscores=True,force_unique_taxon_labels=False) for data, base in zip(datas, bases)]

