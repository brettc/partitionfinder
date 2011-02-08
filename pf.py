import logging
log = logging.getLogger("main")

from optparse import OptionParser
import sys, os

from partfinder import Configuration, ConfigurationError

def load_configuration(pth):
    """Read in the config file and get all the info ready to process. Raise an
    exception if there are any errors"""

    # Get an absolute path so the user is not confused about where we are
    # looking
    pth = os.path.abspath(pth)
    if not os.path.exists(pth) or not os.path.isdir(pth):
        log.error("No such folder: '%s'", pth)
        raise ConfigurationError

    log.info("Using folder: '%s'", pth)
    config = Configuration(pth)
    config.load()
    return config

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
        config = load_configuration(args[0])
        process_configuration(config)
    except ConfigurationError:
        return 1

    # Successful exit
    return 0


if __name__ == "__main__":
    # Well behaved unix programs exits with 0 on success...
    sys.argv = ['arg', '-v', 'example']
    sys.exit(main())

# import logging

# # set up logging to file - see previous section for more details
# logging.basicConfig(level=logging.DEBUG,
                    # format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    # datefmt='%m-%d %H:%M',
                    # filename='/temp/myapp.log',
                    # filemode='w')
# # define a Handler which writes INFO messages or higher to the sys.stderr
# console = logging.StreamHandler()
# console.setLevel(logging.INFO)
# # set a format which is simpler for console use
# formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
# # tell the handler to use this format
# console.setFormatter(formatter)
# # add the handler to the root logger
# logging.getLogger('').addHandler(console)

# # Now, we can log to the root logger, or any other logger. First the root...
# logging.info('Jackdaws love my big sphinx of quartz.')

# # Now, define a couple of other loggers which might represent areas in your
# # application:

# logger1 = logging.getLogger('myapp.area1')
# logger2 = logging.getLogger('myapp.area2')

# logger1.debug('Quick zephyrs blow, vexing daft Jim.')
# logger1.info('How quickly daft jumping zebras vex.')
# logger2.warning('Jail zesty vixen who grabbed pay from quack.')
# logger2.error('The five boxing wizards jump quickly.')
