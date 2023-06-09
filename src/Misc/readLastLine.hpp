#ifndef READLASTLINE_HPP
#define READLASTLINE_HPP

#include <cstdio>

long getFileSize(FILE *fp);
int readLastLine(FILE *fp, char *buffer, int bufferLength);
int readNextIterationNumber(const char *fname);

#endif
