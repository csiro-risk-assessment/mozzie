// Copyright (c) 2024 Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../csvparser.h"

int main(int argc, char **argv)
{
  if (argc != 7)
  {
    printf("You must call this program with the correct arguments.  They are:\n");
    printf(" - conversion_type: Can be either ascii2binary or binary2ascii\n");
    printf(" - nx: The number of grid cells in the x direction\n");
    printf(" - ny: The number of grid cells in the y direction\n");
    printf(" - filetype: Can be inactive_active, wind_raw, wind_processed, generic_float\n");
    printf(" - input file: The ascii-plaintext file to be converted to binary\n");
    printf(" - output file: This will be the binary version of the input file\n");
    return 1;
  }

  int c_from = 0;
  int c_to = 1;
  if (strcmp(argv[1], "ascii2binary") == 0)
  {
    c_from = 0;
    c_to = 1;
  }
  else if (strcmp(argv[1], "binary2ascii") == 0)
  {
    c_from = 1;
    c_to = 0;
  }
  else
  {
    printf("Incorrectly specified conversion_type\n");
  }
    
  
  int err = 0;
  // data, in the same form as in the somewhat obscure cython version.  Hopefully this makes the code clearer!
  char **header = (char**)malloc(sizeof(char*));
  size_t *header_length = (size_t*)malloc(sizeof(size_t));
  unsigned **uint = (unsigned**)malloc(sizeof(unsigned*));
  float **float_data = (float**)malloc(sizeof(float*));
  size_t *wind_length = (size_t*)malloc(sizeof(size_t));

  // read data
  size_t nx = strtoul(argv[2], NULL, 10);
  size_t ny = strtoul(argv[3], NULL, 10);
  if (strcmp(argv[4], "inactive_active") == 0)
  {
    err = parseBool(argv[5], header, header_length, uint, nx, nx * ny, c_from);
  }
  else if (strcmp(argv[4], "wind_raw") == 0)
  {
    err = parseFloat(argv[5], header, header_length, float_data, 2 * nx, 2 * nx * ny, c_from);
  }
  else if (strcmp(argv[4], "wind_processed") == 0)
  {
    err = parseProcessedWind(argv[5], header, header_length, uint, float_data, wind_length, 3, c_from);
  }
  else if (strcmp(argv[4], "generic_float") == 0)
  {
    err = parseFloat(argv[5], header, header_length, float_data, nx, nx * ny, c_from);
  }
  else
  {
    printf("filetype incorrectly specified\n");
    return 1;
  }
  if (err != 0)
  {
    printf("ERROR: parsing returned error code %d\n", err);
    return err;
  }

  // write data
  if (strcmp(argv[4], "inactive_active") == 0)
  {
    err = writeBool(argv[6], header[0], uint[0], nx, nx * ny, c_to);
  }
  else if (strcmp(argv[4], "wind_raw") == 0)
  {
    err = writeFloat(argv[6], header[0], float_data[0], 2 * nx, 2 * nx * ny, c_to);
  }
  else if (strcmp(argv[4], "wind_processed") == 0)
  {
    err = writeProcessedWind(argv[6], header[0], uint[0], float_data[0], wind_length[0], c_to);
  }
  else if (strcmp(argv[4], "generic_float") == 0)
  {
    err = writeFloat(argv[6], header[0], float_data[0], nx, nx * ny, c_to);
  }
  if (err != 0)
  {
    printf("ERROR: writing returned error code %d\n", err);
  }
  return err;
}
