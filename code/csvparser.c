#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "csvparser.h"

int getHeader(FILE *fptr, char **header, size_t *header_length)
{
  // the file is read line-by-line
  // each line is max_line_length long (only this many bytes minus 1 will be read per line)
  size_t max_line_length = 16384;
  char line[max_line_length];
  // Since max_line_length is actually quite large, initialise header[0] to its size
  header[0] = (char*) malloc(max_line_length * sizeof(char));
  if (header[0] == NULL)
  {
    return 1;
  }
  strcpy(header[0], "");
  // max number of chars (without the trailing \0) that can be kept in header[0]
  size_t max_chars_in_header0 = max_line_length - 1;

  while ((fgets(line, sizeof(line), fptr) != NULL) && line[0] == '#') // only process header lines
  {
    header_length[0] = strlen(header[0]);
    while (strlen(line) + header_length[0] > max_chars_in_header0) // the strcat, below, will overflow header[0]
    {
      // allocate a bigger string
      max_chars_in_header0 += max_line_length;
      char * new_header = (char*) malloc((max_chars_in_header0 + 1) * sizeof(char));
      if (new_header == NULL)
      {
        return 1;
      }
      strcpy(new_header, header[0]);
      free(header[0]);
      header[0] = new_header;
    }
    strcat(header[0], line);
  }
  header_length[0] = strlen(header[0]);
  fclose(fptr);
  return 0;
}

int parse(const char *filename, char **header, size_t *header_length)
{
  FILE *fptr = fopen(filename, "r");
  if (fptr == NULL)
  {
    return 1;
  }
  int header_error = getHeader(fptr, header, header_length);
  if (header_error != 0)
  {
    return header_error;
  }
  return 0;
}

