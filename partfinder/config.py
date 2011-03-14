import logging
log = logging.getLogger("config")
from logging.handlers import RotatingFileHandler

import os, shutil

import parser

class ConfigurationError(Exception):
    pass

def _check_file(pth):
    if not os.path.exists(pth) or not os.path.isfile(pth):
        log.error("No such file: '%s'", pth)
        raise ConfigurationError

def _check_folder(pth):
    if not os.path.exists(pth) or not os.path.isdir(pth):
        log.error("No such folder: '%s'", pth)
        raise ConfigurationError

def _make_folder(pth):
    if os.path.exists(pth):
        if not os.path.isdir(pth):
            log.error("Cannot create folder '%s'", pth)
            raise ConfigurationError
    else:
        os.mkdir(pth)

class Settings(object):
    """This holds the user configuration info"""
    def __init__(self):
        pass

    def __getattr__(self, name):
        if name not in self.__dict__:
            log.error("The setting '%s' is not defined. "
                      "Did you initialise the configuration?", 
                      name)
            raise ConfigurationError

        return self.__dict__[name]

class Data(object):
    """Used to hold global data"""
    def __init__(self):
        pass

    # def load_example(self):
        # pass

init_done = False
settings = Settings()
data = Data()

def initialise(pth, force_restart=False):
    global init_done
    if init_done:
        log.error("Cannot initialise more than once")
        raise ConfigurationError

    # Allow for user and environment variables
    pth = os.path.expanduser(pth)
    pth = os.path.expandvars(pth)
    pth = os.path.normpath(pth)
    settings.base_path = pth
    settings.force_restart = force_restart

    _check_folder(settings.base_path)
    log.info("Using folder: '%s'", settings.base_path)

    create_debug_log()

    settings.output_path = os.path.join(settings.base_path, "output")
    if settings.force_restart:
        if os.path.exists(settings.output_path):
            log.warning("Deleting all previous workings in '%s'", settings.output_path)
            shutil.rmtree(settings.output_path)

    _make_folder(settings.output_path)

    # Setup the testing path
    # TODO Should really just run a bunch of tests with --run-tests option
    # How do we do this with nose?
    # settings.test_path = os.path.join(get_root_install_path(), 'tests')

    find_program()

    init_done = True

def get_root_install_path():
    pth = os.path.abspath(__file__)
    # Split off the name and the directory...
    pth, notused = os.path.split(pth)
    pth, notused = os.path.split(pth)
    return pth

def create_debug_log():
    """Add a full debug log file in the folder"""
    settings.log_path = os.path.join(settings.base_path, "partition_finder.log")
    log.info("Full log is in: '%s'", settings.log_path)

    # Append to the log file. we'll get multiple runs then
    log_output = RotatingFileHandler(
        settings.log_path,
        maxBytes=100*1024, # 100K will do?
        backupCount=5,
    )
    log_output.setLevel(logging.DEBUG)
    formatter = logging.Formatter(
        '%(levelname)-8s | %(asctime)s | %(name)-10s | %(message)s',
        datefmt="%Y-%m-%d %H:%M:%S")
    log_output.setFormatter(formatter)
    logging.getLogger('').addHandler(log_output)
    # Mark a new session so it is easy to read in the log file
    log.debug("------------------ NEW LOG SESSION BEGINS -----------------")

def load():
    settings.config_path = os.path.join(settings.base_path, "partition_finder.cfg")
    _check_file(settings.config_path)

    log.info("Loading configuration at '%s'", settings.config_path)
    try:
        p = parser.Parser(settings)
        p.parse_file(settings.config_path)
    except parser.ParserError, p:
        # Catch any parsing errors, print something out, then raise a
        # configuration error, as this is the general error we expect from
        # this part of the process
        log.error(p.format_message())
        raise ConfigurationError

    settings.alignment_path = os.path.join(settings.base_path,
                                        settings.alignment_file)
    _check_file(settings.alignment_path)

def find_program():
    """Locate the binary ..."""

    # TODO: This is a bit crap... maybe look at how other's do it
    # We should really just try and run it. Can we just run it to see what
    # version it is?
    #
    # try:
        # p = Popen(cmd, stdout=PIPE, stderr=STDOUT)
        # out = p.communicate()[0].decode()
        # for k, v in CC_SIGNATURE.items():
            # m = v.search(out)
            # if m:
                # return k
    # except OSError:
        # pass
    # return None

    program_name = 'phyml'
    # if sys.platform == 'win32':
        # program_name += ".exe"

    # Now go back down into programs...
    pth = os.path.join(get_root_install_path(), "programs", program_name)
    pth = os.path.normpath(pth)

    log.debug("Checking for program %s", program_name)
    _check_file(pth)
    log.debug("Found program %s at '%s'", program_name, pth)
    settings.program_path = pth

def remove_tempdir(pth):
    log.debug("Removing temp folder %s", pth)
    shutil.rmtree(pth)

def initialise_temp():
    import tempfile
    import atexit
    tmp = tempfile.mkdtemp()
    atexit.register(remove_tempdir, tmp)
    initialise(tmp, True)

def initialise_example():
    # NOTE: this overwrites everything!
    example_path = os.path.join(get_root_install_path(), 'example')
    initialise(example_path, True)
    load()

def report_settings():
    log.debug("Settings are as follows:")
    for x in settings.__dict__:
        if not x.startswith('__'):
            log.debug("%s: %s", x, getattr(settings, x))

if __name__ == '__main__':
    import logging
    logging.basicConfig(level=logging.DEBUG)
    initialise_temp()
    report_settings()
    
