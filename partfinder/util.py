import logging
log = logging.getLogger("util")
import os


# Base error class
class PartitionFinderError(Exception):
    pass

def check_file_exists(pth):
    if not os.path.exists(pth) or not os.path.isfile(pth):
        log.error("No such file: '%s'", pth)
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
