#ifndef SGC_H__
#define SGC_H__
#ifdef SUPPORT_GZIP_COMPRESSED
#include "zlib.h"
#define fptr gzFile
#define INPUT_OPEN gzopen
#define INPUT_SEEK gzseek
#define INPUT_CLOSE gzclose
#define INPUT_GETS(str, n, stream) gzgets(stream, str, n)
#else
#define fptr FILE *
#define INPUT_OPEN fopen
#define INPUT_SEEK fseek
#define INPUT_CLOSE fclose
#define INPUT_GETS(str, n, stream) fgets(str, n, stream)
#endif
#endif
