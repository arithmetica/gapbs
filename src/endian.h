/* Take from Galois
  Nov 28 2018
 */
#ifndef GALOIS_ENDIAN_H
#define GALOIS_ENDIAN_H

#include <stdint.h>

#define USE_NAIVE_BYTE_SWAP

#ifdef HAVE_ENDIAN_H
# include <endian.h>
#endif

// NB: Wrap these standard functions with different names because
// sometimes le64toh and such are implemented as macros and we don't
// want any nasty surprises.
static inline uint64_t convert_le64(uint64_t x) {
#if !defined(HAVE_BIG_ENDIAN)
  return x;
#elif defined(USE_NAIVE_BYTE_SWAP) || !defined(HAVE_LE64TOH)
  return ((x<<56) & 0xFF00000000000000) | 
         ((x<<40) & 0x00FF000000000000) |
         ((x<<24) & 0x0000FF0000000000) |
         ((x<<8 ) & 0x000000FF00000000) |
         ((x>>8 ) & 0x00000000FF000000) |
         ((x>>24) & 0x0000000000FF0000) |
         ((x>>40) & 0x000000000000FF00) |
         ((x>>56) & 0x00000000000000FF);
#else
  return le64toh(x);
#endif
}

static inline uint32_t convert_le32(uint32_t x) {
#if !defined(HAVE_BIG_ENDIAN)
  return x;
#elif defined(USE_NAIVE_BYTE_SWAP) || !defined(HAVE_LE64TOH)
  return ((x<<24) & 0xFF000000) |
         ((x<<8 ) & 0x00FF0000) |
         ((x>>8 ) & 0x0000FF00) |
         ((x>>24) & 0x000000FF);
#else
  return le32toh(x);
#endif
}

#endif
