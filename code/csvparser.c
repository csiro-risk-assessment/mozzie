#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "csvparser.h"

// we read the entire file in one gulp, to attempt to speed I/O
#define MAX_FILE_LENGTH (100 * 1024 * 1024)
char file_contents[MAX_FILE_LENGTH];

// initial guess at the maximum length of the header.  the header string gets dynamically resized if the file header is longer than this
#define MAX_HEADER_LINE_LENGTH (4 * 1024)

// a single line within the file, got by using strtok_r on file_contents
char * single_line;

// for splitting the file_contents into lines
char * next_line;

int getHeader(char **header, size_t *header_length)
{
  // the file is read line-by-line, using strtok_r
  next_line = file_contents;
  // initialise header[0] to MAX_HEADER_LINE_LENGTH
  header[0] = (char*) malloc(MAX_HEADER_LINE_LENGTH * sizeof(char));
  if (header[0] == NULL)
  {
    return 1;
  }
  strcpy(header[0], "");
  // max number of chars (without the trailing \0) that can be kept in header[0]
  size_t max_chars_in_header0 = MAX_HEADER_LINE_LENGTH - 1;

  while ((single_line = strtok_r(next_line, "\n", &next_line)) && single_line[0] == '#') // only process header lines
  {
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

int read_and_get_header(const char *filename, char **header, size_t *header_length)
{
  FILE *fptr = fopen(filename, "r");
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
  int header_error = getHeader(header, header_length);
  if (header_error != 0)
  {
    return header_error;
  }
  return 0;
}

int getDataBool(unsigned **uint, size_t num_in_row, size_t expected_size)
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
  // go line-by-line through the file_contents
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
      uint[0][num_found] = strtoul(next_data, 0L, 10);
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

  if (num_found != expected_size)
  {
    return 4;
  }
  return 0;
}

int parseBool(const char *filename, char **header, size_t *header_length, unsigned **uint, size_t num_in_row, size_t expected_size)
{
  int oerr = read_and_get_header(filename, header, header_length);
  if (oerr != 0)
  {
    return oerr;
  }
  return getDataBool(uint, num_in_row, expected_size);
}

int parseFloat(const char *filename, char **header, size_t *header_length, float **float_data, size_t num_in_row, size_t expected_size)
{
  int oerr = read_and_get_header(filename, header, header_length);
  if (oerr != 0)
  {
    return oerr;
  }
  return 0;
}

int parseWind(const char *filename, char **header, size_t *header_length, float **wind_data, size_t num_in_row, size_t expected_size)
{
  int oerr = read_and_get_header(filename, header, header_length);
  if (oerr != 0)
  {
    return oerr;
  }
  return 0;
}

int parseProcessedWind(const char *filename, char **header, size_t *header_length, float **processed_wind_data, size_t num_in_row, size_t expected_size)
{
  int oerr = read_and_get_header(filename, header, header_length);
  if (oerr != 0)
  {
    return oerr;
  }
  return 0;
}

