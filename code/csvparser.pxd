# Copyright (c) 2024 Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
cdef extern from "csvparser.c":
    pass

cdef extern from "csvparser.h":
    int parseBool(const char *filename, char **header, size_t *header_length, unsigned **uint, size_t num_in_row, size_t expected_size, int binary)
    int parseFloat(const char *filename, char **header, size_t *header_length, float **float_data, size_t num_in_row, size_t expected_size, int binary)
    int parseProcessedWind(const char *filename, char **header, size_t *header_length, unsigned **uint, float **float_data, size_t *num_found, size_t num_in_row, int binary)
    int writeBool(const char *filename, const char *header, const unsigned *uint, size_t num_in_row, size_t u_length, int binary)
    int writeFloat(const char *filename, const char *header, const float *float_data, size_t num_in_row, size_t f_length, int binary)
    int writeProcessedWind(const char *filename, const char *header, const unsigned *uint, const float *float_data, size_t wind_length, int binary)
