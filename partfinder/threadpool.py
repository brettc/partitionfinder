#
# Modified and vastly simplified from here
# http://code.activestate.com/recipes/203871/
#
import threading
from time import sleep
import sys, os

# Taken from the multiprocessing library
def cpu_count():
    '''
    Returns the number of CPUs in the system
    '''
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
        return num
    # This will have to do
    return 1

_cpus = cpu_count()

class Pool:
    def __init__(self, tasks, numthreads=None):
        """Initialize the thread pool with numthreads workers and all tasks"""
        self.tasks = tasks
        self.tasklock = threading.Condition(threading.Lock())
        self.threads = []
        if numthreads is None:
            numthreads = _cpus
        for i in range(numthreads):
            t = Thread(self)
            self.threads.append(t)
            t.start()

    def next_task(self):
        self.tasklock.acquire()
        try:
            if self.tasks == []:
                return (None, None, None)
            else:
                return self.tasks.pop(0)
        finally:
            self.tasklock.release()
    
    def join(self):
        # Wait for tasks to finish
        while self.tasks != []:
            sleep(.1)

class Thread(threading.Thread):
    def __init__(self, pool):
        threading.Thread.__init__(self)
        self.pool = pool
        
    def run(self):
        while 1:
            cmd, args, callback = self.pool.next_task()
            # If there's nothing to do, just sleep a bit
            if cmd is None:
                break
            elif callback is None:
                cmd(args)
            else:
                callback(cmd(args))
    
# Usage example
if __name__ == "__main__":
    import subprocess, shlex

    def command_task(command):
        print 'starting', command
        p = subprocess.Popen(
            shlex.split(command),
            shell=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)

        # Capture the output, we might put it into the errors
        stdout, stderr = p.communicate()
        print 'done', command
        return stdout

    def callback(output):
        print 'done', len(output)

    # Insert tasks into the queue and let them run
    tasks = [
    (command_task, 'find /Users/brett/Dropbox -name "B*"', callback),
    (command_task, 'find /Users/brett/Dropbox -name "Brett*"', callback),
    (command_task, 'find /Users/brett/Dropbox -name "Brett*"', callback),
    (command_task, 'find /Users/brett/Library -name "Brett*"', callback),
    (command_task, 'find /Users/brett/Dropbox -name "Brett*"', callback),
    (command_task, 'find /Users/brett/Dropbox -name "Brett*"', callback),
    (command_task, 'find /Users/brett/Dropbox -name "Brett*"', callback),
    (command_task, 'find /Users/brett/Dropbox -name "Brett*"', callback),
    ]

    pool = Pool(tasks)
    # print command_task('find /Users/brett/Dropbox -name "Brett*"')
    # When all tasks are finished, allow the threads to terminate
    pool.join()
