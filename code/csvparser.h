#ifndef CSVPARSER_H
#define CSVPARSER_H

int getHeader(FILE *fptr, char **header, size_t *header_length);
int parse(const char *filename, char **header, size_t *header_length);

#endif



