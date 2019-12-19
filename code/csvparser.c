#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "csvparser.h"

int getHeader(const char *filename, char **header, size_t *header_length)
{
  size_t max_line_length = 16384;
  header[0] = (char*) malloc(max_line_length * sizeof(char));
  size_t current_max_header_length = max_line_length - 1;
  char line[max_line_length - 1];
  strcpy(header[0], "");

  FILE *fptr = fopen(filename, "r");
  if (fptr == NULL)
  {
    return 1;
  }

  while ((fgets(line, sizeof(line), fptr) != NULL) && line[0] == '#')
  {
    header_length[0] = strlen(header[0]);
    if (strlen(line) + header_length[0] > current_max_header_length)
    {
      return 1;
    }
    strcat(header[0], line);
  }
  header_length[0] = strlen(header[0]);
  fclose(fptr);
  return 0;
}

