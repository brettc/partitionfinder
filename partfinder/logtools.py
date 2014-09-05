import logging
import os
import re
import textwrap
import inspect

_log_depth = 0
_max_width = 80
_tab_width = 2

# These should be the same size as the _tab_width
_bullet       = ""
_continuation = "..."

def get_logger(fname=None):
    """Pass in the __file__"""

    # Magically get the filename from the calling function
    if fname is None:
        caller_frame = inspect.stack()[1][0]
        fname = caller_frame.f_globals['__file__']

    # Strip the beginning and the extension
    head_tail = os.path.split(fname)
    base_ext = os.path.splitext(head_tail[1])

    log_name = base_ext[0]

    # Trim it max 10 characters
    log_name = log_name[:10]

    # Now wrap it and return it
    # return DumbLogger(logging.getLogger(log_name))
    return SmartLogger(logging.getLogger(log_name))

class DumbLogger(object):
    def __init__(self, logger):
        self.log = logger

    def debug(self, *args):
        self.log.debug(*args)

    def info(self, *args):
        self.log.info(*args)

    def warning(self, *args):
        self.log.warning(*args)

    def error(self, *args):
        self.log.error(*args)

    def push():
        pass

    def pop():
        pass


class SmartLogger(object):
    def __init__(self, logger):
        self.log = logger

    def debug(self, *args):
        msg = self.compose_message(*args)
        self.post_message(msg, self.log.debug)

    def info(self, *args):
        msg = self.compose_message(*args)
        self.post_message(msg, self.log.info)

    def warning(self, *args):
        msg = self.compose_message(*args)
        self.post_message(msg, self.log.warning)

    def error(self, *args):
        msg = self.compose_message(*args)
        self.post_message(msg, self.log.error)

    def format_message(self, msg):
        """Strip multiline comments down to a single line"""
        # First, get rid of tabs and newlines
        warning = re.sub('\s', ' ', msg)

        # Now get rid of all extra spaces
        # http://stackoverflow.com/questions/1546226/
        # the-shortest-way-to-remove-multiple-spaces-in-a-string-in-python
        return ' '.join(warning.split())

    def compose_message(self, *args):
        if len(args) > 1:
            msg = args[0] % args[1:]
        else:
            msg = args[0]
        msg = self.format_message(msg)
        return msg

    def normal_post_message(self, msg, log_function):
        indent_amount = _log_depth * (_tab_width + 1)
        spaces = " " * indent_amount
        log_function(spaces + msg)

    def clever_post_message(self, msg, log_function):
        # How much room have we got?
        indent_amount = _log_depth * (_tab_width + 1)
        local_max_width = _max_width - indent_amount
        spaces = " " * indent_amount

        # Note: if this happens, you're doing too much formatting. Get rid of
        # some of the "push" calls
        # TODO: you can't do an assert here.
        assert local_max_width > 30

        if len(msg) <= local_max_width:
            log_function(spaces + _bullet + msg)
            return

        # Add an extra for the hanging indent
        lines = textwrap.wrap(msg, local_max_width)
        line_iterator = iter(lines)
        first = line_iterator.next()
        log_function(spaces + _bullet + first)
        for next_line in line_iterator:
            log_function(spaces + _continuation + next_line)

    post_message = normal_post_message

    def push(self):
        global _log_depth
        _log_depth += 1

    def pop(self):
        global _log_depth
        _log_depth -= 1


class indented(object):
    def __init__(self, logger=None, msg=None):
        self.logger = logger
        self.msg = msg

    def __enter__(self):
        if self.logger and self.msg:
            self.logger.info(self.msg)
        self.logger.push()

    def __exit__(self, type, value, traceback):
        self.logger.pop()


class log_info(object):
    """Decorator for wrapping functions indented"""
    def __init__(self, logger, msg):
        self.logger = logger
        self.msg = msg

    def __call__(self, fn):
        def indented_fn(*args, **kwargs):
            with indented(self.logger, self.msg):
                fn(*args, **kwargs)
        return indented_fn
