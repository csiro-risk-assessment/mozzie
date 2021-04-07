#ifndef CSVPARSER_H
#define CSVPARSER_H

// Private method.  ASSUMING that read_and_get_header has been called, in order to file the file_contents string, place the header string into header, and return header_length too.  After calling this function, file_contents will have NULL characters in place of \n for all header lines (those at the start of the file that start with '#') and next_line will point to the first byte after the header data.
int getHeader(char **header, size_t *header_length);

// Private method.  Open file, read all its data into file_contents, close the file, and then call getHeader to populate "header" and "header_length".  After calling this function, file_contents will have NULL characters in place of \n for all header lines (those at the start of the file that start with '#') and next_line will point to the first byte after the header data.
int read_and_get_header(const char *filename, char **header, size_t *header_length, int binary);

// Private method.  ASSUMING getHeader has just been called, so that next_line is now pointing at the data inside the file_contents string, extract the remaining data (assumed to be unsigned longs - either 0 or 1) into uint.
int getDataBool(unsigned **uint, size_t header_len, size_t num_in_row, size_t expected_size, int binary);

// Private method.  ASSUMING getHeader has just been called, so that next_line is now pointing at the data inside the file_contents string, extract the remaining data (assumed to be floats) into float_data.
int getDataFloat(float **float_data, size_t header_len, size_t num_in_row, size_t expected_size, int binary);

// Private method.  ASSUMING getHeader has just been called, so that next_line is now pointing at the data inside the file_contents string, extract the remaining data (assumed to be unsigned, unsigned, float) into uint and float_data, putting the number of data lines found in num_found[0]
int getProcessedWind(unsigned **uint, size_t header_len, float **float_data, size_t *num_found, size_t num_in_row, int binary);

// Public method.  Parse data in filename, placing header data (as a single string) into "header", and boolean (unsigned, either 0 or 1) into uint.  If binary=0 the file is assumed ascii CSV; if binary=1 the file is a binary representation of this, created by ab_convert (very little error-checking is done if binary=1).  Returns 0 upon no error.  Error codes are:
// 1: cannot open or close file, or allocate header string (memory error)
// 2: length of file exceeds MAX_FILE_LENGTH
// 3: header is too long in the file, and cannot allocate memory for a new, longer, header string.  Or, cannot allocate uint
// 4: number of data does not equal expected_size
// 5: each row does not contain num_in_row data
// 6: data is not 0 or not 1
int parseBool(const char *filename, char **header, size_t *header_length, unsigned **uint, size_t num_in_row, size_t expected_size, int binary);

// Public method.  Parse data in filename, placing header data (as a single string) into "header", and float data into float_data.  If binary=0 the file is assumed ascii CSV; if binary=1 the file is a binary representation of this, created by ab_convert (very little error-checking is done if binary=1).  Returns 0 upon no error.  Error codes are:
// 1: cannot open or close file, or allocate header string (memory error)
// 2: length of file exceeds MAX_FILE_LENGTH
// 3: header is too long in the file, and cannot allocate memory for a new, longer, header string.  Or, cannot allocate float_data
// 4: number of data does not equal expected_size
// 5: each row does not contain num_in_row data
int parseFloat(const char *filename, char **header, size_t *header_length, float **float_data, size_t num_in_row, size_t expected_size, int binary);

// Public method.  Parse data in filename, placing header data (as a single string) into "header", and cell numbers into uint and probabilities into float_data.  The number of data rows found is put into num_found[0].   If binary=0 the file is assumed ascii CSV; if binary=1 the file is a binary representation of this, created by ab_convert (very little error-checking is done if binary=1).  Returns 0 upon no error.  Error codes are:
// 1: cannot open or close file, or allocate header string (memory error)
// 2: length of file exceeds MAX_FILE_LENGTH
// 3: header is too long in the file, and cannot allocate memory for a new, longer, header string.  Or, cannot allocate float_data
// 4: number of data does not equal expected_size
// 5: each row does not contain num_in_row data
int parseProcessedWind(const char *filename, char **header, size_t *header_length, unsigned **uint, float **float_data, size_t *num_found, size_t num_in_row, int binary);

// PublicMethod.  Write the active/inactive information in uint to filename, with given header.  num_in_row = number of grid cells in the x direction.  u_length = length of uint.  binary controls whether the output is binary or ascii-CSV
int writeBool(const char *filename, const char *header, const unsigned *uint, size_t num_in_row, size_t u_length, int binary);

// PublicMethod.  Write the floats to filename, with given header.  num_in_row = number of grid cells in the x direction.  f_length = length of float_data.  binary controls whether the output is binary or ascii-CSV
int writeFloat(const char *filename, const char *header, const float *float_data, size_t num_in_row, size_t f_length, int binary);

// PublicMethod.  Write the uints and floats to filename, with given header.  wind_length is the number of entries in float_data (there must be double this number in uint).  binary controls whether the output is binary or ascii-CSV
int writeProcessedWind(const char *filename, const char *header, const unsigned *uint, const float *float_data, size_t wind_length, int binary);

#endif



