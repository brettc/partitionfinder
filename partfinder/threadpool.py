#
# Modified and vastly simplified from here
# http://code.activestate.com/recipes/203871/
#
#Copyright (C) 2011 Robert Lanfear and Brett Calcott
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
#program and the PyParsing library both of which are protected by their
#own licenses and conditions, using PartitionFinder implies that you
#agree with those licences and conditions as well.

import logging
log = logging.getLogger("threadpool")
import threading
from time import sleep
import sys, os

# Catch these exceptions
from util import PhylogenyProgramError

# Taken from the multiprocessing library
def cpu_count():
    if sys.platform == 'win32':
        try:
            num = int(os.environ['NUMBER_OF_PROCESSORS'])
        except (ValueError, KeyError):
            num = 0
    elif sys.platform == 'darwin':
        try:
            num = int(os.popen('sysctl -n hw.ncpu').read())
        except ValueError:
            num = 0
    else:
        try:
            num = os.sysconf('SC_NPROCESSORS_ONLN')
        except (ValueError, OSError, AttributeError):
            num = 0

    if num > 1:
        log.info("You appear to have %s cpus", num)
        return num
    # This will have to do
    log.info("I cannot detect any more than one processor...")
    return 1

_cpus = cpu_count()

class Pool(object):
    def __init__(self, tasks, numthreads=-1):
        """Initialize the thread pool with numthreads workers and all tasks"""
        self.more_tasks = True
        self.tasks = tasks
        self.task_lock = threading.Condition(threading.Lock())
        self.threads = []

        numtasks = len(tasks)
        if numtasks == 0:
            log.warning("You did not give any tasks to do...")
            self.more_tasks = False
            return

        if numthreads <= 1:
            numthreads = _cpus
        if numtasks < numthreads:
            numthreads = numtasks

        log.debug("Creating %s threads for %s tasks", numthreads, numtasks)
        for i in range(numthreads):
            t = Thread(self)
            self.threads.append(t)
            t.start()

    def next_task(self):
        self.task_lock.acquire()
        try:
            if self.tasks == []:
                self.more_tasks = False
                return None, None
            else:
                return self.tasks.pop(0)
        finally:
            self.task_lock.release()

    def join(self):
        # TODO: I don't think we need this bit....
        # Wait till all tasks have been taken
        while self.more_tasks:
            sleep(.1)
        # ... now wait for them all to finish
        for t in self.threads:
            t.join()

class Thread(threading.Thread):
    def __init__(self, pool):
        threading.Thread.__init__(self)
        self.pool = pool

    def run(self):
        while 1:
            cmd, args = self.pool.next_task()
            # If there's nothing to do, return
            if cmd is None:
                break
            try:
                cmd(*args)
            except PhymlError:
                # Catch the ones we know about, the error should already have
                # been reported. Stop operation though.
                break


