#
# Modified and vastly simplified from here
# http://code.activestate.com/recipes/203871/
#
import logging
log = logging.getLogger("threadpool")
import threading
from time import sleep
import sys, os

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

    if num >= 1:
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
            result = cmd(*args)

