// Copyright 2013 Daniel Kang

#include <stdint.h>

#ifndef CONSTS_H_
#define CONSTS_H_

#ifndef NUM_IT
  #define NUM_IT (12)
#endif

#ifndef KBIG
  #define KBIG (8)
#endif
#ifndef K
  #define K (1000)
#endif
#ifndef K_BETA
  #define K_BETA (200)
#endif
#ifndef NUM_COV
  #define NUM_COV (1)
#endif
#ifndef RESOL
  #define RESOL (4)
#endif

#define K2 (2*K)
#define K2_BETA (2 * K_BETA)
#define NUM_BETA (NUM_COV * K2_BETA)
#define BUF_PAD (K2 * 10)

#define ALLOC_SIZE(SIZE) ((SIZE) + BUF_PAD)
#define KALLOC (K2 / RESOL)

#define YM (1 << (KBIG*2))
#define YGRAD_SIZE (KALLOC * YM)
#define YALL_SIZE ((KALLOC * ((1 << (2*KBIG + 2)) - 1)) / 3)

#define DEFINE_CONSTS(KS) \
  const int __attribute__((unused)) M = (1 << (KS*2));\
  const int __attribute__((unused)) XSIZE_GRAD = (KALLOC*M);\
  const int __attribute__((unused)) XSIZE_ALL = ((KALLOC * ((1 << (2*KS + 2)) - 1)) / 3);

// #define SAFE_THREADING

typedef float   ftype;
typedef double  dtype;
typedef int32_t gtype;

#endif  // CONSTS_H_
