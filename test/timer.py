from __future__ import print_function

import time

class Timer(object):
    """
    Time context for benchmarking
    """
    def __init__(self, name=""):
        self.name = name
        self._start = None
        self._end = None

    def start(self):
        self._start = time.time()
        self._end = self._start

    def stop(self):
        self._end = time.time()

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, *args):
        self.stop()
        if self.name:
            print("%s : %0.3f seconds" % (self.name, self.elapsed))

    @property
    def elapsed(self):
        if self._end is None:
            end_time = time.time()
        else:
            end_time = self._end
        return end_time - self._start

    def __str__(self):
        return "Timer(elapsed=%s)" % self.elapsed

def benchmark(f, n_repeats=3, warmup=True, name=""):
    """
    Run the given function f repeatedly, return the average elapsed time.
    """
    if warmup:
        f()

    total_time = 0
    for i in range(n_repeats):
        iter_name = "%s (iter #%d)" % (name, i + 1,)
        with Timer(iter_name) as t:
            f()
        total_time += t.elapsed
    return total_time / n_repeats
