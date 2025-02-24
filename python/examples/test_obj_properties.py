class Agent:
    def __init__(self):
        self.propa = None
        self.propb = None

    def calc(self, input):
        out = input*self.propa
        self.propb = out
        return out
    
    def calcdirect(self, input, propa):
        return input*propa
    
if __name__ == "__main__":
    import time
    import cProfile, pstats, io
    from pstats import SortKey

    a = Agent()
    a.propa = 5
    nmax = 1000000 

    tic = time.time()
    for k in range(nmax):
        a.calc(3)
    toc = time.time()
    dur = toc - tic

    tic = time.time()
    val = a.propa
    for k in range(nmax):
        a.calcdirect(3, val)
    toc = time.time()
    dur2 = toc - tic

    print(f"property: {dur:.3f} s, direct: {dur2:.3f} s, ratio: {dur/dur2}")

    # Now repeat with the cProfiler
    pr = cProfile.Profile()
    pr.enable()
    for k in range(nmax):
        a.calc(3)
    pr.disable()
    s = io.StringIO()
    sortby = SortKey.TIME
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats()
    print(s.getvalue())
    print(dir(ps))
    print(dir(s))
