#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "csvparser.h"

// we read the entire file in one gulp, to attempt to speed I/O
#define MAX_FILE_LENGTH (100 * 1024 * 1024)
char file_contents[MAX_FILE_LENGTH];

// these are relevant for the header lines only
#define MAX_HEADER_LINE_LENGTH (4 * 1024)
char header_line[MAX_HEADER_LINE_LENGTH];

int getHeader(FILE *fptr, char **header, size_t *header_length)
{
  // the file is read line-by-line
  // each line is MAX_HEADER_LINE_LENGTH long (only this many bytes minus 1 will be read per line)
  // Since MAX_HEADER_LINE_LENGTH is actually quite large, initialise header[0] to its size
  header[0] = (char*) malloc(MAX_HEADER_LINE_LENGTH * sizeof(char));
  if (header[0] == NULL)
  {
    return 1;
  }
  strcpy(header[0], "");
  // max number of chars (without the trailing \0) that can be kept in header[0]
  size_t max_chars_in_header0 = MAX_HEADER_LINE_LENGTH - 1;

  while ((fgets(header_line, sizeof(header_line), fptr) != NULL) && header_line[0] == '#') // only process header lines
  {
    header_length[0] = strlen(header[0]);
    while (strlen(header_line) + header_length[0] > max_chars_in_header0) // the strcat, below, will overflow header[0]
    {
      // allocate a bigger string
      max_chars_in_header0 += MAX_HEADER_LINE_LENGTH;
      char * new_header = (char*) malloc((max_chars_in_header0 + 1) * sizeof(char));
      if (new_header == NULL)
      {
        return 1;
      }
      strcpy(new_header, header[0]);
      free(header[0]);
      header[0] = new_header;
    }
    strcat(header[0], header_line);
  }
  header_length[0] = strlen(header[0]);
  return 0;
}

int parse(const char *filename, char **header, size_t *header_length, unsigned **uint)
{
  FILE *fptr = fopen(filename, "r");
  if (fptr == NULL)
  {
    return 1;
  }
  size_t num_bytes_read = fread(file_contents, 1, MAX_FILE_LENGTH, fptr);
  if (num_bytes_read >= MAX_FILE_LENGTH)
  {
    return 2;
  }
  if (fclose(fptr) != 0)
  {
    return 1;
  }
  fptr = fopen(filename, "r");
  if (fptr == NULL)
  {
    return 1;
  }
  int header_error = getHeader(fptr, header, header_length);
  if (header_error != 0)
  {
    return header_error;
  }
  fclose(fptr);
  return 0;
}

