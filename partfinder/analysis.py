import multiprocessing
import subprocess

import scheme, subset, partition

results_cache = {}

def process(s):
    pass

def multiprocess_all(s):
    pass

def process_all(s):
    pass

def work(cmd):
    # return subprocess.call(cmd, shell=False)
    p = subprocess.Popen('ctags -R .'.split(),
                            shell=False,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)

    stdout, stderr = p.communicate()
    return stdout

def multi():
    count = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=count)
    x = pool.map(work, ['ls'] * 20)
    print x

class Subset:
    # From Subset
    def analyse(self):
        # Check first to see if we've got the results, otherwise calculate and
        # cache them.
        if self.partitions in results_cache:
            log.debug("Returning cached result for %s", self)
            return results_cache[self.partitions]

        log.debug("Calculating result for %s", self)
        result = self._really_analyse()
        results_cache[self.partitions] = result
        return result

    def _really_analyse(self):
        fname = self.make_filename()
        sa = alignment.SubsetAlignment(
            fname, config.settings.source_alignment, self)
        return sa.analyse()

class X:
    def source_exists(self):
        return os.path.exists(self.source_path)
    def analysis_exists(self):
        return os.path.exists(self.analysis_path)

    def same_as_saved(self):
        spec, slen = self.read_source()
        if spec == self.species:
            log.debug("Are same")
            return True
        else:
            log.warning("different")
            return False

    def analyse(self):
        if self.source_exists():
            # We're already written a file of this name
            log.debug("%s already exists at '%s'", self, self.source_path)
            if self.same_as_saved():
                log.debug("%s Same", self)
                same_as = True
                XXXXXXXX
        else:
            # Otherwise write it out
            self.write_source()

        fresh_analysis = True
        if self.analysis_exists():
            log.debug("Reading in previous analysis of %s", self)
            output = file(self.analysis_path, 'r').read()
            fresh_analysis = False
        else:
            output = modelgen.run(self.source_path)
            log.debug("Saving analysis output of %s to %s", self,
                      self.analysis_path)
            open(self.analysis_path, 'w').write(output)

        log.debug("Parsing ModelGenerator output for %s", self)
        result = modelgen.parse(output)

        if fresh_analysis:
            # Show them how long it took.
            log.debug("New analysis of %s took %d seconds", self, result.processing_time)
        return result
    @property
    def source_path(self):
        return self.path + "." + alignment_format

    @property
    def analysis_path(self):
        return self.path + ".out"

    @property
    def path(self):
        # Defer to something that get's defined differently in subclasses
        # And cache it.
        if not hasattr(self, "_path"):
            self._path = self.get_path()
        return self._path

    def get_path(self):
        # Don't use this class -- use one of the subclasses below
        raise NotImplemented

if __name__ == '__main__':

"""
# 1. Try multiprocessing using pool

# 2. QUEUES
#
def worker():
    while True:
        item = q.get()
        do_work(item)
        q.task_done()

q = Queue()
for i in range(num_worker_threads):
     t = Thread(target=worker)
     t.daemon = True
     t.start()

for item in source():
    q.put(item)

q.join()       # block until all tasks are done

# 3. threading + subprocess
import threading
import subprocess

class MyClass(threading.Thread):
    def __init__(self):
        self.stdout = None
        self.stderr = None
        threading.Thread.__init__(self)

    def run(self):
        p = subprocess.Popen('rsync -av /etc/passwd /tmp'.split(),
                             shell=False,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)

        self.stdout, self.stderr = p.communicate()

myclass = MyClass()
myclass.start()
myclass.join()
print myclass.stdout
"""
