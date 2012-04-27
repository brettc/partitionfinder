
# alternative curry-by-closure, Alex Martelli, 2001/06/28
# Nick Perkins' latest version of curry-by-closure (from c.l.p, for the
# record) is more general (no named args at all, thus no accidental name
# capture -- except, I think, for args and create_time_kwds...) and quite
# readable, although marginally more verbose, due to good name choices for
# intermediate locals. It also has curry-from-left (for positional args) and
# call-time-dominates (for named args) semantics, which may be popular: 

def curry(*args, **create_time_kwds):
    func = args[0]
    create_time_args = args[1:]
    def curried_function(*call_time_args, **call_time_kwds):
        args = create_time_args + call_time_args
        kwds = create_time_kwds.copy()
        kwds.update(call_time_kwds)
        return func(*args, **kwds)
    return curried_function

# half-finished currying statement, which allows for partial replacement
# of positional args
def named_curry(f, **create_time_kwds):
    func = f
    argpos = {}
    kwargs = {}
    for k, v in create_time_kwds.items():
        if k.startswith('arg'):
            argnum = k[3:]
            pos = int(argnum)
            argpos[pos] = v
        else:
            kwargs[k] = v
    
    # FIXME need to check the args here

    # create_time_args = args[1:]
    def curried_function(*call_time_args, **call_time_kwds):
        # now here, need to interleave ...
        pos = 1
        args = []
        call_time_args = list(call_time_args[::-1]) # reversed, so we can pop
        while 1:
            posarg = argpos.get(pos, None)
            # should pop?
            if posarg:
                args.append(posarg)
            else:
                # run out?
                if not call_time_args:
                    break
                args.append(call_time_args.pop())
            pos += 1

        # FIXME, check that we have used all create_time_args
        kwds = kwargs.copy()
        kwds.update(call_time_kwds)
        return func(*args, **kwds)

    return curried_function


def X(a, b, c, blarg=11):
    print a, b, c, blarg


ac = named_curry(X, arg1=2, arg2=3, blarg=1)
ac(4)
