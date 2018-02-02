/*
 * Written 1998 by Andreas Johansson <ajo@wopr.campus.luth.se>
 *
 * Use this software at your own risk, the author of it will not take
 * any responsibility for any damage caused by it.
 */

#include <math.h>
#include "pow_nint.h"

double pow_nint_tab[PN_SIZE];

void init_pow_nint (void)
{
  int i;
  
  for (i = 1; i < PN_SIZE; i++) {
    pow_nint_tab[i] = pow ((double)i - 0.4054, 4.0/3.0);
  }
}
