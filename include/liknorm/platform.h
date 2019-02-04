#ifndef LIKNORM_PLATFORM_H
#define LIKNORM_PLATFORM_H

#include <assert.h>

#ifdef _WIN32
#ifdef LIKNORM_API_EXPORTS
#define LIKNORM_API __declspec(dllexport)
#else
#define LIKNORM_API __declspec(dllimport)
#endif
#else
#define LIKNORM_API
#endif

/* Borrowed from GLIB. */
#if __GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1)
#define LIKNORM_DEPR __attribute__((__deprecated__))
#else
#define LIKNORM_DEPR
#endif

#ifndef static_assert
#define static_assert(...)
#endif

#endif /* LIKNORM_PLATFORM_H */
