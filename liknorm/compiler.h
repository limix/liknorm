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

#endif
