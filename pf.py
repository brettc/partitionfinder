import logging
log = logging.getLogger("main")

from optparse import OptionParser
import sys, os

from partfinder import Configuration, ConfigurationError

def main():
    usage = """usage: python %prog <folder_containing_config>

    README.txt should be here
    """
    parser = OptionParser(usage)
    # parser.add_option("-f", "--file", dest="filename",
                      # help="read data from FILENAME")
    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="Show verbose (debug) output")

    options, args = parser.parse_args()
    # We should have one argument -- the folder to read the configuration from
    if not args:
        parser.print_help()
        return 2

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
        # Any exceptions and we fail
        return 1

    # Successful exit
    return 0

if __name__ == "__main__":
    # Well behaved unix programs exits with 0 on success...
    sys.argv = ['arg', '-v', 'example']
    sys.exit(main())


