/*******************************************************************************
    PRODIGAL (PROkaryotic DynamIc Programming Genefinding ALgorithm)
    Copyright (C) 2007-2016 University of Tennessee / UT-Battelle

    Code Author:  Doug Hyatt

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*******************************************************************************/

#ifndef BITMAP_H_
#define BITMAP_H_
#include <stddef.h>

/* Test a bit, 0 = not set, 1 = set */
static inline unsigned char test(unsigned char *bm, size_t ndx) {
  return ( bm[ndx>>3] & (1 << (ndx&0x07))?1:0 ); 
}

/* Clear a bit (set it to 0) */
static inline void clear(unsigned char *bm, size_t ndx) {
  bm[ndx>>3] &= ~(1 << (ndx&0x07));
}

/* Set a bit to 1 */
static inline void set(unsigned char *bm, size_t ndx) {
  bm[ndx>>3] |= (1 << (ndx&0x07));
}

/* Flip a bit's value 0->1 or 1->0 */
static inline void toggle(unsigned char *bm, size_t ndx) {
  bm[ndx>>3] ^= (1 << (ndx&0x07)); 
}

/* Get two bits */
static inline unsigned char nuc(unsigned char *bm, size_t ndx) {
  ndx *= 2;
  return ((bm[ndx >> 3] >> (ndx & 0x07)) & 0x03);
}

/* Get six bits, which can be in two bytes
   out-of-bound read may happen but will not affect result */
static inline unsigned char trinuc(unsigned char *bm, size_t ndx) {
  ndx *= 2;
  unsigned byte = ndx >> 3;
  unsigned stretch = (bm[byte + 1] << 8) | bm[byte];
  return ((stretch >> (ndx & 0x07)) & 0x3F);
}

/* Get twelve bits, which can be in four bytes
   mer_ndx(6) is relatively common and should be replaced with this
   (shorter assembly).  unfortunately the order is different, which
   breaks tinf->gene_dc compatibility! */
static inline unsigned char hexnuc(unsigned char *bm, int ndx) {
  ndx *= 2;
  unsigned byte = ndx >> 3;
  unsigned stretch = (bm[byte + 3] << 24) | (bm[byte + 2] << 16) |
                     (bm[byte + 1] << 8) | bm[byte];
  return ((stretch >> (ndx & 0x07)) & 0x3F);
}
#endif
