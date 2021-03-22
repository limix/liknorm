#include "report.h"
#include <stdio.h>

void liknorm_error(char const *err, ...)
{
    va_list params;
    va_start(params, err);
    fprintf(stderr, err, params);
    fputc('\n', stderr);
    va_end(params);
}
