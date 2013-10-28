cdef extern:
    void setup_c()

cdef extern:
    void discus_loop()

def setup():
    setup_c()

def interactive():
    discus_loop()
