import numpy as	np

cdef extern int diff_power

cdef extern:
    void setup_c()

cdef extern:
    void discus_loop_c()

def setup():
    setup_c()

def test():
    print diff_power

def discus_loop():
    discus_loop_c()
