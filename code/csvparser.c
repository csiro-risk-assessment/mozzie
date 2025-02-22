// Copyright (c) 2024 Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "csvparser.h"

#if defined(_WIN32) || defined(_WIN64)
/* We are on Windows */
# define strtok_r strtok_s
#endif

// we read the entire file in one gulp, to attempt to speed I/O
#define MAX_FILE_LENGTH (200 * 1024 * 1024) // may have to be increased to read large files, or decreased because of computer-memory restrictions
char file_contents[MAX_FILE_LENGTH];

// initial guess at the maximum length of the header.  the header string gets dynamically resized if the file header is longer than this
#define MAX_HEADER_LINE_LENGTH (4 * 1024)

// a single line within the file, got by using strtok_r on file_contents
char * single_line;

// for splitting the file_contents into lines
char * next_line;

// additional header information written into binary files
char binary_header[] = "#strlen of remaining header = ";

int getHeader(char **header, size_t *header_length, int binary)
{
  // the file is read line-by-line, using strtok_r

  size_t max_chars_in_header0 = MAX_HEADER_LINE_LENGTH - 1;   // max number of chars (without the trailing \0) that can be kept in header[0]

  size_t *actual_header_strlen = (size_t*) malloc(sizeof(size_t));
  if (binary == 0)
  {
    // initialise header[0] to MAX_HEADER_LINE_LENGTH
    header[0] = (char*) malloc(MAX_HEADER_LINE_LENGTH * sizeof(char));
    max_chars_in_header0 = MAX_HEADER_LINE_LENGTH - 1;
    // point next_line at the start of the file
    next_line = file_contents;
  }
  else
  {
    memcpy(actual_header_strlen, file_contents + strlen(binary_header) * sizeof(char), sizeof(size_t));
    // initialise header[0] to the length just read from the file
    header[0] = (char*) malloc(((*actual_header_strlen) + 1) * sizeof(char));
    max_chars_in_header0 = *actual_header_strlen;
    // point next_line at the start of the header
    next_line = file_contents + strlen(binary_header) * sizeof(char) + sizeof(size_t);
  }

  if (header[0] == NULL)
  {
    return 1;
  }
  strcpy(header[0], "");

  single_line = next_line;
  while (*next_line == '#' && single_line != NULL && (binary == 0 || strlen(header[0]) < (*actual_header_strlen)))
  {
    single_line = strtok_r(next_line, "\n", &next_line);
    header_length[0] = strlen(header[0]);
    while (strlen(single_line) + 1 + header_length[0] > max_chars_in_header0) // the strcats, below, will overflow header[0]
    {
      // allocate a bigger string
      max_chars_in_header0 += MAX_HEADER_LINE_LENGTH;
      char * new_header = (char*) malloc((max_chars_in_header0 + 1) * sizeof(char));
      if (new_header == NULL)
      {
        return 3;
      }
      strcpy(new_header, header[0]);
      free(header[0]);
      header[0] = new_header;
    }
    strcat(header[0], single_line);
    strcat(header[0], "\n");
  }
  header_length[0] = strlen(header[0]);
  return 0;
}

int read_and_get_header(const char *filename, char **header, size_t *header_length, int binary)
{
  FILE *fptr = NULL;
  if (binary == 0)
    fptr = fopen(filename, "r");
  else if (binary == 1)
    fptr = fopen(filename, "rb");
  if (fptr == NULL)
  {
    return 1;
  }
  size_t num_bytes_read = fread(file_contents, 1, MAX_FILE_LENGTH, fptr);
  if (num_bytes_read >= MAX_FILE_LENGTH - 1) // -1 because of the \0 character added below
  {
    return 2;
  }
  if (fclose(fptr) != 0)
  {
    return 1;
  }
  file_contents[num_bytes_read] = '\0';
  int header_error = getHeader(header, header_length, binary);
  if (header_error != 0)
  {
    return header_error;
  }
  return 0;
}

int getDataBool(unsigned **uint, size_t header_len, size_t num_in_row, size_t expected_size, int binary)
{
  uint[0] = (unsigned*) malloc(expected_size * sizeof(unsigned));
  if (uint[0] == NULL)
  {
    return 3;
  }

  // number of data found in file_contents
  size_t num_found = 0;
  // number of data found in single_line
  size_t num_found_in_line = 0;

  char * next_data;
  if (binary == 0)
  {
    // go line-by-line through the file_contents
    single_line = strtok_r(next_line, "\n", &next_line);
    do
    {
      if (!single_line)
      {
	continue;
      }
      num_found_in_line = 0;
      while((next_data = strtok_r(single_line, ",", &single_line)))
      {
	if (num_found >= expected_size)
        {
	  return 4;
        }
        uint[0][num_found] = strtoul(next_data, NULL, 10);
        if (!(uint[0][num_found] == 0 || uint[0][num_found] == 1))
        {
	  return 6;
        }
        num_found += 1;
        num_found_in_line += 1;
      }
      if (num_found_in_line != num_in_row)
      {
        return 5;
      }
    } while ((single_line = strtok_r(next_line, "\n", &next_line)));
  }
  else
  {
    memcpy(uint[0], file_contents + strlen(binary_header) * sizeof(char) + sizeof(size_t) + header_len * sizeof(char), expected_size * sizeof(unsigned));
    num_found = expected_size;
  }

  if (num_found != expected_size)
  {
    return 4;
  }
  return 0;
}

int getDataFloat(float **float_data, size_t header_len, size_t num_in_row, size_t expected_size, int binary)
{
  float_data[0] = (float*) malloc(expected_size * sizeof(float));
  if (float_data[0] == NULL)
  {
    return 3;
  }

  // number of data found in file_contents
  size_t num_found = 0;
  // number of data found in single_line
  size_t num_found_in_line = 0;

  char * next_data;
  if (binary == 0)
  {
    // go line-by-line through the file_contents
    single_line = strtok_r(next_line, "\n", &next_line);
    do
    {
      if (!single_line)
      {
        continue;
      }
      num_found_in_line = 0;
      while((next_data = strtok_r(single_line, ",", &single_line)))
      {
        if (num_found >= expected_size)
        {
	  return 4;
        }
        float_data[0][num_found] = strtof(next_data, NULL);
        num_found += 1;
        num_found_in_line += 1;
      }
      if (num_found_in_line != num_in_row)
      {
	return 5;
      }
    } while ((single_line = strtok_r(next_line, "\n", &next_line)));
  }
  else
  {
    memcpy(float_data[0], file_contents + strlen(binary_header) * sizeof(char) + sizeof(size_t) + header_len * sizeof(char), expected_size * sizeof(float));
    num_found = expected_size;
  }

  if (num_found != expected_size)
  {
    return 4;
  }
  return 0;
}

int getProcessedWind(unsigned **uint, size_t header_len, float **float_data, size_t *num_found, size_t num_in_row, int binary)
{
  uint[0] = (unsigned*) malloc(2 * MAX_FILE_LENGTH * sizeof(unsigned));
  float_data[0] = (float*) malloc(MAX_FILE_LENGTH * sizeof(float));
  if (uint[0] == NULL || float_data[0] == NULL)
  {
    return 3;
  }

  // number of lines found in file_contents
  num_found[0] = 0;
  // number of data found in single_line
  size_t num_found_in_line = 0;

  char * next_data;
  if (binary == 0)
  {
    // go line-by-line through the file_contents
    single_line = strtok_r(next_line, "\n", &next_line);
    do
    {
      if (!single_line)
      {
	continue;
      }
      num_found_in_line = 0;
      while((next_data = strtok_r(single_line, ",", &single_line)))
      {
        if (num_found[0] >= MAX_FILE_LENGTH)
        {
	  return 2;
        }
        num_found_in_line += 1;
        if (num_found_in_line == 1)
        {
          uint[0][2 * num_found[0]] = strtoul(next_data, NULL, 10);
        }
        else if (num_found_in_line == 2)
        {
          uint[0][2 * num_found[0] + 1] = strtoul(next_data, NULL, 10);
        }
        else if (num_found_in_line == 3)
        {
          float_data[0][num_found[0]] = strtof(next_data, NULL);
        }
      }
      num_found[0] += 1;
      if (num_found_in_line != num_in_row)
      {
        return 5;
      }
    } while ((single_line = strtok_r(next_line, "\n", &next_line)));
  }
  else
  {
    char * startat = file_contents + strlen(binary_header) * sizeof(char) + sizeof(size_t) + header_len * sizeof(char);
    memcpy(num_found, startat, sizeof(size_t));
    memcpy(uint[0], startat + sizeof(size_t), 2 * num_found[0] * sizeof(unsigned));
    memcpy(float_data[0], startat + sizeof(size_t) + 2 * num_found[0] * sizeof(unsigned), num_found[0] * sizeof(float));
  }

  return 0;
}

int parseBool(const char *filename, char **header, size_t *header_length, unsigned **uint, size_t num_in_row, size_t expected_size, int binary)
{
  int oerr = read_and_get_header(filename, header, header_length, binary);
  if (oerr != 0)
  {
    return oerr;
  }
  return getDataBool(uint, header_length[0], num_in_row, expected_size, binary);
}

int parseFloat(const char *filename, char **header, size_t *header_length, float **float_data, size_t num_in_row, size_t expected_size, int binary)
{
  int oerr = read_and_get_header(filename, header, header_length, binary);
  if (oerr != 0)
  {
    return oerr;
  }
  return getDataFloat(float_data, header_length[0], num_in_row, expected_size, binary);
}

int parseProcessedWind(const char *filename, char **header, size_t *header_length, unsigned **uint, float **float_data, size_t *num_found, size_t num_in_row, int binary)
{
  int oerr = read_and_get_header(filename, header, header_length, binary);
  if (oerr != 0)
  {
    return oerr;
  }
  return getProcessedWind(uint, header_length[0], float_data, num_found, num_in_row, binary);
}

void writeBinaryHeader(FILE *fptr, const char *header)
{
    fwrite(binary_header, sizeof(char), strlen(binary_header), fptr);
    size_t len_header = strlen(header);
    size_t * len_ptr = &len_header;
    fwrite(len_ptr, sizeof(size_t), 1, fptr);
}

int writeBool(const char *filename, const char *header, const unsigned *uint, size_t num_in_row, size_t u_length, int binary)
{
  size_t i, row, ind;
  size_t num_rows = u_length / num_in_row;
  FILE *fptr = NULL;
  if (binary == 0)
    fptr = fopen(filename, "w");
  else if (binary == 1)
    fptr = fopen(filename, "wb");
  if (fptr == NULL)
  {
    return 1;
  }
  if (binary == 0)
  {
    fprintf(fptr, "%s", header);
    ind = 0;
    for (row = 0; row < num_rows; ++row)
    {
      for (i = 0; i < num_in_row - 1; ++i)
      {
	fprintf(fptr, "%u,", uint[ind]);
	ind += 1;
      }
      fprintf(fptr, "%u\n", uint[ind]);
      ind += 1;
    }
  }
  else if (binary == 1)
  {
    writeBinaryHeader(fptr, header);
    fwrite(header, sizeof(char), strlen(header), fptr);
    fwrite(uint, sizeof(unsigned), u_length, fptr);
  }
  if (fclose(fptr) != 0)
  {
    return 1;
  }
  return 0;
}

int writeFloat(const char *filename, const char *header, const float *float_data, size_t num_in_row, size_t f_length, int binary)
{
  size_t i, row, ind;
  size_t num_rows = f_length / num_in_row;
  FILE *fptr = NULL;
  if (binary == 0)
    fptr = fopen(filename, "w");
  else if (binary == 1)
    fptr = fopen(filename, "wb");
  if (fptr == NULL)
  {
    return 1;
  }
  if (binary == 0)
  {
    fprintf(fptr, "%s", header);
    ind = 0;
    for (row = 0; row < num_rows; ++row)
    {
      for (i = 0; i < num_in_row - 1; ++i)
      {
	fprintf(fptr, "%f,", float_data[ind]);
	ind += 1;
      }
      fprintf(fptr, "%f\n", float_data[ind]);
      ind += 1;
    }
  }
  else if (binary == 1)
  {
    writeBinaryHeader(fptr, header);
    fwrite(header, sizeof(char), strlen(header), fptr);
    fwrite(float_data, sizeof(float), f_length, fptr);
  }
  if (fclose(fptr) != 0)
  {
    return 1;
  }
  return 0;
}


int writeProcessedWind(const char *filename, const char *header, const unsigned *uint, const float *float_data, size_t wind_length, int binary)
{
  size_t row;
  FILE *fptr = NULL;
  if (binary == 0)
    fptr = fopen(filename, "w");
  else if (binary == 1)
    fptr = fopen(filename, "wb");
  if (fptr == NULL)
  {
    return 1;
  }
  if (binary == 0)
  {
    fprintf(fptr, "%s", header);
    for (row = 0; row < wind_length; ++row)
    {
      fprintf(fptr, "%u,%u,%f\n", uint[2 * row], uint[2 * row + 1], float_data[row]);
    }
  }
  else if (binary == 1)
  {
    writeBinaryHeader(fptr, header);
    fwrite(header, sizeof(char), strlen(header), fptr);
    fwrite(&wind_length, sizeof(size_t), 1, fptr);
    fwrite(uint, sizeof(unsigned), 2 * wind_length, fptr);
    fwrite(float_data, sizeof(float), wind_length, fptr);
  }
  if (fclose(fptr) != 0)
  {
    return 1;
  }
  return 0;
}
