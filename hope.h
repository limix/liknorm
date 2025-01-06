#ifndef HOPE_H
#define HOPE_H

#include "support.h"

#define CLOSE(actual, desired)                                                 \
    CLOSE2((double)(actual), (double)(desired), __HOPE_REL_TOL(actual), 0.0)

#define CLOSE2(actual, desired, rel_tol, abs_tol)                              \
    __hope_close2((double)(actual), (double)(desired), (double)(rel_tol),      \
                  (double)(abs_tol), __FILE__, __LINE__)

#define EQ(actual, desired) __hope_eq((actual), (desired), __FILE__, __LINE__)

#define ISNULL(a) __hope_isnull((a) != NULL, #a, __FILE__, __LINE__)

#define NOTNULL(a) __hope_notnull((a) == NULL, #a, __FILE__, __LINE__)

#define COND(a) __hope_cond(#a, !(a), __FILE__, __LINE__)

inline static int hope_status(void) { return __hope_errors; }

#endif
