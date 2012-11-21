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
#program, the RAxML program, the PyParsing library, and the python-cluster library 
#all of which are protected by their own licenses and conditions, using 
#PartitionFinder implies that you agree with those licences and conditions as well.

import logging
log = logging.getLogger("util")
import os
import fnmatch


# Base error class
class PartitionFinderError(Exception):
    pass


class PhylogenyProgramError(PartitionFinderError):
    pass


def check_file_exists(pth):
    if not os.path.exists(pth) or not os.path.isfile(pth):
        if pth.count("partition_finder.cfg") > 0:
            log.error("Failed to find configuration file: '%s'. "
                      "For PartitionFinder to run, there must be a file called 'partition_finder.cfg' "
                      "located in the same folder as your alignment. Please check and try again.", pth)
            raise PartitionFinderError
        else:
            log.error(
                "Failed to find file: '%s'. Please check and try again.", pth)
            raise PartitionFinderError


def check_folder_exists(pth):
    if not os.path.exists(pth) or not os.path.isdir(pth):
        log.error("No such folder: '%s'", pth)
        raise PartitionFinderError


def get_root_install_path():
    pth = os.path.abspath(__file__)
    # Split off the name and the directory...
    pth, not_used = os.path.split(pth)
    pth, not_used = os.path.split(pth)
    return pth


def make_dir(pth):
    if os.path.exists(pth):
        if not os.path.isdir(pth):
            log.error("Cannot create folder '%s'", pth)
            raise PartitionFinderError
    else:
        os.mkdir(pth)


def remove_runID_files(aln_pth):
    """remove all files that match a particular run_ID. Useful for cleaning out directories
    but ONLY after a whole analysis of a subset is completely finished, be careful!"""
    head, tail = os.path.split(aln_pth)
    run_ID = os.path.splitext(tail)[0]
    head = os.path.abspath(head)
    fnames = os.listdir(head)
    fs = fnmatch.filter(fnames, '*%s*' % run_ID)
    for f in fs:
        try:
            os.remove(os.path.join(head, f))
        except OSError:
            # Don't complain if you can't delete them
            # (This is here because we sometimes try and delete things twice in
            # the threading).
            pass
