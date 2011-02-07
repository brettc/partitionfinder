import logging
log = logging.getLogger("config")

from optparse import OptionParser
import sys

def main():
    usage = """usage: python %prog <folder_containing_config>

    README.txt should be here
    """
    parser = OptionParser(usage)
    # parser.add_option("-f", "--file", dest="filename",
                      # help="read data from FILENAME")
    parser.add_option("-d", "--debug",
                      action="store_true", dest="debug",
                      help="show debug output")

    options, args = parser.parse_args()
    if options.debug:
        level = logging.DEBUG
    else:
        level = logging.INFO


    # Now evaluate

    logging.basicConfig(level=level)
    log.debug("Now showing all debug messages...")

if __name__ == "__main__":
    main()
