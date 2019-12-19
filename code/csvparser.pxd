cdef extern from "csvparser.c":
    pass

cdef extern from "csvparser.h":
    int getHeader(const char *filename, char *header)
