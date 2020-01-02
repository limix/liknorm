#ifndef REPORT_H
#define REPORT_H

#include "liknorm.h"
#include <stdarg.h>

#if defined(HAVE_ATTR_FORMAT)
#define ATTR_FORMAT __attribute__((format(printf, 1, 2)))
#else
#define ATTR_FORMAT
#endif

void liknorm_error(char const *err, ...) ATTR_FORMAT;

#endif
