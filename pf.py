import logging
log = logging.getLogger("main")

from optparse import OptionParser
import sys, os

from partfinder import Configuration, ConfigurationError, ProcessingError

def main():
    usage = """usage: python %prog [-vc] <foldername>

    README.txt should be here
    """
    parser = OptionParser(usage)
    parser.add_option(
        "-v", "--verbose",
        action="store_true", dest="verbose",
        help="show verbose (debug) output")
    parser.add_option(
        "-c", "--check-only",
        action="store_true", dest="check_only",
        help="Just check the configuration files, don't do any processing")
    # Other options...?
    # --force-restart
    #
    options, args = parser.parse_args()

    # We should have one argument: the folder to read the configuration from
    if not args:
        # Otherwise exit, printing the help
        parser.print_help()
        return 2

    if options.verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO

    logging.basicConfig(
        format='%(levelname)-8s | %(message)s',
        level=level
    )

    # Load, using the first argument as the folder
    try:
        config = Configuration(args[0])
        config.load()
    except ConfigurationError:
        log.error("Configuration Failure: Please correct problems and rerun")
        # Any exceptions and we fail
        return 1

    log.info("Configuration appears to be okay.")
    if options.check_only:
        log.info("Exiting without processing as requested...")
        return 0

    # Now try processing everything....
    log.info("Beginning processing.")
    try:
        config.process()
    except ProcessingError:
        log.error("Processing Error!")
        return 1

    # Successful exit
    log.info("Success: processing complete.")
    return 0

if __name__ == "__main__":
    # Well behaved unix programs exits with 0 on success...
    # sys.argv = ['arg', '-v', 'example']
    # sys.argv = ['arg', 'example']
    # sys.argv = ['arg', '-vc', 'example']
    sys.exit(main())


