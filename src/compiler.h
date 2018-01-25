#ifndef COMPILER_H
#define COMPILER_H

#define LIKELY(x) x
#define UNLIKELY(x) x

// Branch prediction hints
#if defined(__has_builtin)
#if __has_builtin(__builtin_expect)
#undef LIKELY
#define LIKELY(x) __builtin_expect(x, 1)
#undef UNLIKELY
#define UNLIKELY(x) __builtin_expect(x, 0)
#endif
#endif

#if defined(_MSC_VER)
#if _MSC_VER <= 1500
#include "hcephes.h"
#include <float.h>
#define log1p hcephes_log1p
inline double fmax(double left, double right) {
  return (left > right) ? left : right;
}

inline double fmin(double left, double right) {
  return (left < right) ? left : right;
}

#define lgamma hcephes_lgam
#define isnan(x) _isnan(x)
#define isfinite(x) _finite(x)
#endif
#endif

#endif
