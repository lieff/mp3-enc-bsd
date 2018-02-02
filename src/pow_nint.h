/*
 * Written 1998 by Andreas Johansson <ajo@wopr.campus.luth.se>
 *
 * Use this software at your own risk, the author of it will not take
 * any responsibility for any damage caused by it.
 */

#include <math.h>

#define PN_SIZE 4096

extern double pow_nint_tab[PN_SIZE];
extern void init_pow_nint (void);

static inline int pow_nint (double x)
{
  int step, pos, p;

  pos = 1; p = 0; step = 1;
  
  while (pos < (PN_SIZE>>1) ) {
    if (x < pow_nint_tab[pos])
      break;
    
    p = pos;
    pos += step;
    step <<= 1;
  }
  step >>= 1;
  pos -= step;
  step >>= 1;
  
  if (step) {
    while (step) {
      if (x < pow_nint_tab[pos])
	pos -= step;
      else {
	p = pos;
	pos += step;
      }
      
      step >>= 1;
    }
    if (x >= pow_nint_tab[pos])
      p = pos;
  }
  
  return (p);
}

