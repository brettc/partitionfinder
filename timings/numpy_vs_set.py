import timeit
setup_fmt = """
import random
from hashlib import md5
from cPickle import dumps
from numpy import union1d, concatenate, array

cols = range({0})
x = random.sample(cols, {1})
y = random.sample(cols, {1})
x1 = set(x)
y1 = set(y)
x2 = array(x)
y2 = array(y)

def meth0():
    q = x1 | y1
    c = list(q)
    c.sort()
    x = md5(dumps(c, -1)).hexdigest()
    return q, c, x

def meth1():
    q = union1d(x2, y2)
    # q = concatenate((x2, y2))
    # q.sort()
    x = md5(q.tostring()).hexdigest()
    return q, x
"""

if __name__ == '__main__':
    rep = 3
    number = 100
    for column_size in 100, 1000, 10000, 100000:
        print 'using array size\t{0} ----'.format(column_size)
        sample_size = column_size/10
        setup = setup_fmt.format(column_size, sample_size)
        t_set = min(timeit.repeat('meth0()', setup=setup, repeat=rep, number=number))
        t_numpy = min(timeit.repeat('meth1()', setup=setup, repeat=rep, number=number))
        print "using sets\t{:.6f}\nusing numpy\t{:.6f}".format(t_set, t_numpy)


