import logging
log = logging.getLogger("config")

import os
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

def get_root_install_path():
    pth = os.path.abspath(__file__)
    # Split off the name and the directory...
    pth, not_used = os.path.split(pth)
    pth, not_used = os.path.split(pth)
    return pth

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

settings = Settings()

def load(pth):
    # Allow for user and environment variables
    pth = os.path.expanduser(pth)
    pth = os.path.expandvars(pth)
    pth = os.path.normpath(pth)
    # pth = os.path.abspath(pth)

    settings.base_path = pth
    _check_folder(settings.base_path)
    log.info("Using folder: '%s'", settings.base_path)

    settings.config_path = os.path.join(
        settings.base_path, "partition_finder.cfg")
    _check_file(settings.config_path)

    settings.analysis_path = os.path.join(settings.base_path, "analysis")

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
                                        settings.alignment)
    _check_file(settings.alignment_path)
    report_settings()

def report_settings():
    log.debug("Settings are as follows:")
    for x in settings.__dict__:
        if not x.startswith('__'):
            log.debug("%s: %s", x, getattr(settings, x))

if __name__ == '__main__':
    import logging
    logging.basicConfig(level=logging.DEBUG)
    report_settings()
    
