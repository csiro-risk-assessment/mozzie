#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "csvparser.h"

int getHeader(const char *filename, char *header)
{
  // probably better to have char**header, and do malloc on header (see http://docs.cython.org/en/latest/src/tutorial/strings.html)
  strcpy(header, "");
  FILE *fptr = fopen(filename, "r");
  if (fptr == NULL)
  {
    return 1;
  }
  char line[1023];
  fgets(line, 1022, (FILE*)fptr);
  if (line[0] == '#')
  {
    strcpy(header, line);
  }
  fclose(fptr);
  return 0;
}

