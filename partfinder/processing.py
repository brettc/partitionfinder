import multiprocessing
import subprocess

def work(cmd):
    # return subprocess.call(cmd, shell=False)
    p = subprocess.Popen('ctags -R .'.split(),
                            shell=False,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)

    stdout, stderr = p.communicate()
    return stdout

if __name__ == '__main__':
    count = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=count)
    x = pool.map(work, ['ls'] * 20)
    print x

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
