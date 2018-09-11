#pragma once

#if defined(__GNUC__)
#   define BENCH_HOT_SPOT [[gnu::hot]]
#   define BENCH_STRONG_INLINE [[gnu::always_inline]] inline
#   define BENCH_WEAK_INLINE inline
#   define BENCH_DONT_INLINE [[gnu::noinline]] inline
#elif defined(__clang__)
#   define BENCH_HOT_SPOT [[gnu::hot]]
#   define BENCH_STRONG_INLINE [[gnu::always_inline]] inline
#   define BENCH_WEAK_INLINE inline
#   define BENCH_DONT_INLINE [[gnu::noinline]] inline
#else
#   define BENCH_HOT_SPOT
#   define BENCH_STRONG_INLINE inline
#   define BENCH_WEAK_INLINE inline
#   define BENCH_DONT_INLINE inline
#endif
