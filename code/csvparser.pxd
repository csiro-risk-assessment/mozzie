cdef extern from "csvparser.c":
    pass

cdef extern from "csvparser.h":
    int parseBool(const char *filename, char **header, size_t *header_length, unsigned **uint, size_t num_in_row, size_t expected_size)
    int parseFloat(const char *filename, char **header, size_t *header_length, float **float_data, size_t num_in_row, size_t expected_size)
    int parseWind(const char *filename, char **header, size_t *header_length, float **wind_data, size_t num_in_row, size_t expected_size)
    int parseProcessedWind(const char *filename, char **header, size_t *header_length, float **processed_wind_data, size_t num_in_row, size_t expected_size)
