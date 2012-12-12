#Copyright (C) 2012 Robert Lanfear and Brett Calcott
#
#This program is free software: you can redistribute it and/or modify it
#under the terms of the GNU General Public License as published by the
#Free Software Foundation, either version 3 of the License, or (at your
#option) any later version.
#
#This program is distributed in the hope that it will be useful, but
#WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#General Public License for more details. You should have received a copy
#of the GNU General Public License along with this program.  If not, see
#<http://www.gnu.org/licenses/>. PartitionFinder also includes the PhyML
#program, the RAxML program, the PyParsing library, and the python-cluster library
#all of which are protected by their own licenses and conditions, using
#PartitionFinder implies that you agree with those licences and conditions as well.

import logging
log = logging.getLogger("config")
import config

# TODO: Not currently used
# Activation should chdir, and maybe do some other stuff
# So maybe need an 'activate' function on the config?
# Should also clear out subsets, in the cache?

class Current(object):
    """Keep a bunch of stuff current, that can be reinitialised"""
    def __init__(self):
        self._config = None

    def activate_config(self, c):
        assert isinstance(c, config.Configuration)

        if self._config is not None:
            log.debug("Resetting old configuration...")
            self._config.reset()

        log.debug("Assigning a new configuration...")
        self._config = c

    @property
    def active_config(self):
        if self._config is None:
            log.error("No configuration is currently active...")

        return self._config


current = Current()
