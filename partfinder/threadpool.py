# Copyright (C) 2012 Robert Lanfear and Brett Calcott
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details. You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# PartitionFinder also includes the PhyML program, the RAxML program, and the
# PyParsing library, all of which are protected by their own licenses and
# conditions, using PartitionFinder implies that you agree with those licences
# and conditions as well.
#
# This code modified and vastly simplified from here
# http://code.activestate.com/recipes/203871/

import logtools
log = logtools.get_logger()

import threading
from time import sleep
import multiprocessing

_cpus = None

def get_cpu_count():
    global _cpus
    if _cpus is not None:
        return _cpus

    try:
        _cpus = multiprocessing.cpu_count()
    except NotImplementedError:
        _cpus = 1
        log.warning('I cannot detect the number of processors...')
    log.info("Using %s cpus", _cpus)

    return _cpus


class Pool(object):
    def __init__(self, tasks, numthreads=-1):
        """Initialize the thread pool with numthreads workers and all tasks"""
        self.more_tasks = True
        self.tasks = tasks
        self.task_lock = threading.Condition(threading.Lock())
        self.threads = []
        self.failed = False

        numtasks = len(tasks)
        if numtasks == 0:
            log.warning("You did not give any tasks to do...")
            self.more_tasks = False
            return

        if numthreads <= 1:
            numthreads = get_cpu_count()
        if numtasks < numthreads:
            numthreads = numtasks

        self.numtasks = numtasks
        self.curtask = 0

        log.debug("Creating %s threads for %s tasks", numthreads, numtasks)
        for i in range(numthreads):
            t = Thread(self)
            self.threads.append(t)
            t.start()

    def next_task(self):
        self.task_lock.acquire()
        try:
            if self.curtask == self.numtasks:
                self.more_tasks = False
                return None, None
            else:
                task = self.tasks[self.curtask]
                self.curtask += 1 
                return task
        finally:
            self.task_lock.release()

    def kill(self, e):
        self.task_lock.acquire()
        self.tasks = []
        self.more_tasks = False
        self.failed = True
        self.exception = e
        self.task_lock.release()

    def join(self):
        # TODO: I don't think we need this bit....
        # Wait till all tasks have been taken
        while self.more_tasks:
            sleep(.001)
        # ... now wait for them all to finish
        for t in self.threads:
            t.join()

        if self.failed:
            raise self.exception


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
            except Exception as e:
                # The error should already have been reported.
                # Stop operation and kill the entire pool. Then reraise the
                # error
                self.pool.kill(e)
                break
