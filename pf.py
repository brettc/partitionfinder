import logging
log = logging.getLogger("main")

from optparse import OptionParser
import sys, os

from partfinder import Configuration, ConfigurationError

def process_configuration(config):
    """Now process all the information we have"""
    log.info("Beginning processing...")

def main():
    usage = """usage: python %prog <folder_containing_config>

    README.txt should be here
    """
    parser = OptionParser(usage)
    # parser.add_option("-f", "--file", dest="filename",
                      # help="read data from FILENAME")
    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="show verbose output")

    options, args = parser.parse_args()
    if not args:
        parser.print_help()
        return 1

    if options.verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO

    logging.basicConfig(level=level)

    # Right now try
    try:
        config = Configuration(args[0])
        config.load()
        process_configuration(config)
    except ConfigurationError:
        # log.error("Configuration error occurred. Please check log")
        return 1


    # Successful exit
    return 0


if __name__ == "__main__":
    # Well behaved unix programs exits with 0 on success...
    sys.argv = ['arg', '-v', 'example']
    sys.exit(main())


