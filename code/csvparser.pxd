cdef extern from "csvparser.c":
    pass

cdef extern from "csvparser.h":
    int parse(const char *filename, char **header, size_t *header_length)
