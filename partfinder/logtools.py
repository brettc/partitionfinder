import logging
import os
import re
import textwrap

_log_depth = 0
_max_width = 80
_tab_width = 4


def get_logger(fname):
    """Pass in the __file__"""

    # Strip the beginning and the extension
    head_tail = os.path.split(fname)
    base_ext = os.path.splitext(head_tail[1])

    log_name = base_ext[0]

    # Now wrap it and return it
    return Logger(logging.getLogger(log_name))


class Logger(object):
    def __init__(self, logger):
        self.log = logger

    def debug(self, msg):
        self.post_message(msg, self.log.debug)

    def info(self, msg):
        self.post_message(msg, self.log.info)

    def warning(self, msg):
        self.post_message(msg, self.log.warning)

    def error(self, msg):
        self.post_message(msg, self.log.error)

    def post_message(self, msg, log_function):
        msg = self.format_message(msg)
        self.piecemeal_message(msg, log_function)

    def piecemeal_message(self, msg, log_function):
        # How much room have we got?
        indent_amount = _log_depth * _tab_width
        local_max_width = _max_width - indent_amount
        spaces = " " * indent_amount

        # Note: if this happens, you're doing too much formatting. Get rid of
        # some of the "push" calls
        assert local_max_width > 30

        if len(msg) <= local_max_width:
            log_function(spaces + msg)
            return

        lines = textwrap.wrap(msg, local_max_width)
        for l in lines:
            log_function(spaces + l)

    def format_message(self, msg):
        """Strip multiline comments down to a single line"""
        # First, get rid of tabs and newlines
        warning = re.sub('\s', ' ', msg)

        # Now get rid of all extra spaces
        # http://stackoverflow.com/questions/1546226/
        # the-shortest-way-to-remove-multiple-spaces-in-a-string-in-python
        return ' '.join(warning.split())

    @staticmethod
    def push():
        global _log_depth
        _log_depth += 1

    @staticmethod
    def pop():
        global _log_depth
        _log_depth -= 1
