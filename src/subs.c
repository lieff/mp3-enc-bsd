/* Replacement FFT routine for mp3 encoding.
   Wedged in by Mike Cheng  :  http://www.cryogen.com/mikecheng
   Based upon split-radix fft by Malvar
       from book:  "Signal Processing with Lapped Transforms". Malvar, HS.
   
   The reason this fft is so damn fast, is because mp3 encoding only uses
   *real* FFTs.  The previous fft routine included in the iso source is for
   a complex->complex transform.  But if the imaginary part of the signal
   is zero (as it is for encoding sound), then you only have to do a 
   real->complex fft.  This makes this new fft routine about twice the
   speed of the old one.

   There is another real->complex fft that is about another 30% than the
   split-radix (in theory), but it may take me a lot longer to put in.

   later
   mike
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define BLKSIZE 1024
#define MAXLOGM  25
#define TWOPI  6.28318530717958647692
#define SQHALF   0.707106781186547524401

void BR_permute();
void srrec();
void srfft();
void rsfft();
void srifft();
void enphi();
void enphinew();
void fft();

void fft(float x_real[BLKSIZE],float x_imag[BLKSIZE], float energy[BLKSIZE], float phi[BLKSIZE], int N)
{
int logm;

if (N==1024)
  logm=10;
else if (N==256)
  logm=8;
 else {
   printf("logm error\n");
   exit(1);
 }
   rsfft(x_real,logm);
   enphinew(x_real,energy,phi,N);
}
void enphinew(float x_real[BLKSIZE],float energy[BLKSIZE], float phi[BLKSIZE], int N)
{
  int i;
  float *ep, *pp, *epn, *ppn;
  float *xp1, *xpn;

  ep = energy;
  pp = phi;
  xp1 = x_real;
  xpn = &x_real[N-1];

   *ep++ = *xp1 * *xp1;
   *pp++ = atan2( 0.0, (double)(*xp1++));

   if (N==1024) {
   for (i=1;i<512;i++) {
     *ep = *xp1 * *xp1 + *xpn * *xpn;
     if (*ep < 0.0005) {
       *ep++ = 0.0005;
       *pp++ = 0.0;
       xpn--;
       xp1++;
     }
     else {
       ep++;
       *pp++ = atan2( -(double)(*xpn--), (double)(*xp1++) );
     }

   }

   ep = &energy[513];
   pp = &phi[513];
   epn = &energy[511];
   ppn = &phi[511];
   for (i=1;i<512;i++) {
     *ep++ = *epn--;
     *pp++ = - *ppn--;
   }

   energy[512] = x_real[512] * x_real[512];
   phi[512] = atan2 (0.0, (double)x_real[512]);
   }
   else {
   for (i=1;i<N/2;i++) {
     *ep = *xp1 * *xp1 + *xpn * *xpn;
     if (*ep < 0.0005) {
       *ep++ = 0.0005;
       *pp++ = 0.0;
       xpn--;
       xp1++;
     }
     else {
       ep++;
       *pp++ = atan2( -(double)(*xpn--), (double)(*xp1++) );
     }

   }

   ep = &energy[N/2+1];
   pp = &phi[N/2+1];
   epn = &energy[N/2-1];
   ppn = &phi[N/2-1];
   for (i=1;i<N/2;i++) {
     *ep++ = *epn--;
     *pp++ = - *ppn--;
   }

   energy[N/2] = x_real[N/2] * x_real[N/2];
   phi[N/2] = atan2 (0.0, (double)x_real[N/2]);
   }
}
/**********************************************************
* Data unshuffling according to bit-reversed indexing.
*
* Bit reversal is done using Evans' algorithm. (Ref: DMW Evans,
* "An Improved Digit-Reversal Permutation Algorithm",
*  IEEE Trans. ASSP, Aug 1987, pp. 1120-25.
*
***********************************************************/

static int brseed[256]; /* Evans' seed table */
static int brsflg;   /* flag for table building */

void BR_permute(float *x, int logm)
{
register int lg2, n;
int i,j,imax;
int off,fj, gno, *brp;
float tmp, *xp, *xq;

lg2 = logm >> 1;
n = 1 << lg2;

if (logm & 1) lg2++;

/* create seed table if not yet built */
if (brsflg != logm) {
  brsflg = logm;
  brseed[0] = 0;
  brseed[1] = 1;
  for(j=2; j <= lg2; j++){
    imax =  1 << (j - 1);
    for (i=0; i < imax; i++){
      brseed[i] <<= 1;
      brseed[i + imax] = brseed[i] + 1;
   }
 }
}

/* unshuffling loop */
for(off = 1; off < n; off++){
fj = n*brseed[off]; i = off; j = fj;
tmp = x[i]; x[i] = x[j]; x[j] = tmp;
xp = &x[i];
brp = &brseed[1];
for(gno = 1; gno < brseed[off]; gno++){
  xp += n;
  j = fj + *brp++;
  xq = x + j ;
  tmp = *xp; *xp = *xq; *xq = tmp;
}
}


}

/**************************************************************
*                                                             *
* Recursive part of split radix FFT algorithm
*
***************************************************************/

void srrec(float *xr, float *xi, int logm)
{

static int m, m2, m4, m8, nel, n;
static float *xr1, *xr2, *xi1, *xi2;
static float *cn, *spcn, *smcn, *c3n, *spc3n, *smc3n;
static float tmp1, tmp2, ang, c, s;
static float *tab[MAXLOGM];

/* check range of logm */

if ((logm < 0) || (logm > MAXLOGM)){
 fprintf(stderr, "Error: SRFFT logm = %d is out of bounds [%d %d]\n", logm, 0, MAXLOGM);
  exit(1);
}

/* compute trivial cases */
if (logm < 3){
if (logm == 2){
xr2 = xr + 2;
xi2= xi + 2;
tmp1 = *xr + *xr2;
*xr2 = *xr - *xr2;
*xr = tmp1;
tmp1 = *xi + *xi2;
*xi2 = *xi - *xi2;
*xi = tmp1;
xr1 = xr + 1;
xi1 = xi + 1;
xr2++;
xi2++;
tmp1 = *xr1 + *xr2;
*xr2 = *xr1 - *xr2;
*xr1 = tmp1;
tmp1 = *xi1 +  *xi2;
*xi2 = *xi1 - *xi2;
*xi1 = tmp1;
xr2 = xr + 1;
xi2 = xi + 1;
tmp1 = *xr + *xr2;
*xr2 = *xr - *xr2;
*xr = tmp1;
tmp1 = *xi + *xi2;
*xi2 = *xi - *xi2;
*xi = tmp1;
xr1 = xr +2;
xi1 = xi+2;
xr2 = xr+3;
xi2 = xi + 3;
tmp1 = *xr1 + *xi2;
tmp2 = *xi1 + *xr2;
*xi1 = *xi1 - *xr2;
*xr2 = *xr1 - *xi2;
*xr1 = tmp1;
*xi2 = tmp2;
return;
}
else if (logm == 1){
xr2 = xr + 1;
xi2 = xi + 1;
tmp1 = *xr + *xr2;
*xr2 = *xr - *xr2;
*xr = tmp1;
tmp1 = *xi + *xi2;
*xi2 = *xi - *xi2;
*xi = tmp1;
return;
}
else if (logm == 0) return;
}

/* compute a few constants */

m = 1 << logm; m2 = m / 2; m4 = m2 / 2; m8 = m4 / 2;


/* build tables of butterfly coefficients if necessary */
if ((logm >=4) && (tab[logm-4] == NULL)){
/* allocate memory for tables */
nel = m4 - 2;
if ((tab[logm-4] = (float *) calloc(6 * nel, sizeof(float)))
    == NULL){
  exit(1);
}



/* initialize pointers */
cn = tab[logm-4]; spcn = cn + nel; smcn = spcn + nel;
c3n = smcn + nel; spc3n = c3n + nel; smc3n = spc3n + nel;


/* compute tables */
for (n = 1; n < m4; n++){
  if (n == m8) continue;
  ang = n * TWOPI / m;
  c = cos(ang); s = sin(ang);
  *cn++ = c; *spcn++ = -(s + c); *smcn++ = s - c;
  ang = 3*n*TWOPI/m;
  c = cos(ang); s = sin(ang);
  *c3n++ = c; *spc3n++ = -(s + c); *smc3n++ = s - c;
}
}


/* step 1 */
xr1 = xr; xr2 = xr1 + m2;
xi1 = xi; xi2 = xi1 + m2;

for (n =0; n < m2; n++){
  tmp1 = *xr1 + *xr2;
  *xr2 = *xr1 - *xr2;
  *xr1 = tmp1;
  tmp2 = *xi1 + *xi2;
  *xi2 = *xi1 - *xi2;
  *xi1 = tmp2;
  xr1++; xr2++; xi1++; xi2++;
} 

/* Step 2 */
xr1 = xr + m2;  xr2 = xr1 + m4;
xi1 = xi + m2;  xi2 = xi1+ m4;
for (n = 0; n < m4; n++){
  tmp1 = *xr1 + *xi2;
  tmp2 = *xi1 + *xr2;
  *xi1 = *xi1 - *xr2;
  *xr2 = *xr1 - *xi2;
  *xr1 = tmp1;
  *xi2 = tmp2;
  xr1++; xr2++; xi1++; xi2++;
}

/* Steps 3&4 */
xr1 = xr + m2; xr2 = xr1 + m4;
xi1 = xi + m2; xi2 = xi1 + m4;

if (logm >= 4) {
  nel = m4 -2;
  cn = tab[logm-4]; spcn = cn + nel; smcn = spcn + nel;
  c3n = smcn + nel; spc3n = c3n + nel; smc3n = spc3n + nel;
}

xr1++; xr2++; xi1++; xi2++;

for(n=1; n < m4; n++){
   if (n == m8){
     tmp1 = SQHALF*(*xr1 + *xi1);
     *xi1 = SQHALF*(*xi1 - *xr1);
     *xr1 = tmp1;
     tmp2 = SQHALF*(*xi2 - *xr2);
     *xi2 = -SQHALF*(*xr2 + *xi2);
     *xr2 = tmp2;
   }
   else {
     tmp2 = *cn++ *(*xr1 + *xi1);
     tmp1 = *spcn++ * *xr1 + tmp2;
     *xr1 = *smcn++ * *xi1 + tmp2;
     *xi1 = tmp1;
     tmp2 = *c3n++ * (*xr2 + *xi2);
     tmp1 = *spc3n++ * *xr2 + tmp2;
     *xr2 = *smc3n++ * *xi2 + tmp2;
     *xi2 = tmp1;
   }
xr1++; xr2++; xi1++; xi2++;
 }


/* call ssrec again with half DFT length */
srrec(xr, xi, logm -1);
/* call ssrec again twice with one quarter DFT length.
Constants have to be recomputed because they are static! */

m = 1 << logm; m2 = m/2;
srrec(xr+m2, xi+m2, logm - 2);

m = 1 << logm; m4 = 3*(m/4);
srrec(xr+m4, xi+m4, logm - 2);
}

/*****************************************************************
*  Direct transform                                               *
*
*******************************************************************/
        
void srfft(float *xr, float *xi, int logm)
{
/* call recursive routine */
srrec(xr, xi, logm);

/* output array unshuffling using bit-reversed indices */
if (logm > 1){
  BR_permute(xr, logm);
  BR_permute(xi, logm);
}
}
/* --------------------------------------------------------------- */
void enphi(float *xrp, float *xip, float *ep, float *pp,int N)
{
  /* put energy / phi calcs in separate routine */
int i;

 for (i=0;i<N;i++) {
    *ep = *xrp * *xrp + *xip * *xip;
    if (*ep <= 0.0005) {
        *pp++ = 0; 
        *ep++ = 0.0005;
        xrp++;
        xip++;
        }
    else { 
        *pp++ = atan2((double)(*xip++), (double)(*xrp++));
        ep++;
        }
    }
}
/* ------------------------------------------------------  */
/* spli8t radix fast fourert transform for real valued inputs.
   malvar's source code */


/* error exit */
static void error_exit(void)
{
  exit(1);
}

/* recursive part of rsfft algo. not externally callable */
static void rsrec(float *x, int logm) {
  static int m,m2,m4,m8,nel, n;
  static float *xr1, *xr2, *xi1;
  static float *cn, *spcn, *smcn;
  static float tmp1, tmp2, ang, c, s;
  static float *tab[MAXLOGM];

  /* check rang of logm */
  if ((logm <0)|| (logm>MAXLOGM)) {
    fprintf(stderr,"log out of range: %i\n",logm);
    error_exit();
  }

  /* compute trivial cases */
  if (logm <2) {
    if (logm ==1) { /* length m = 2*/
      xr2 = x + 1;
      tmp1 = *x + *xr2;
      *xr2 = *x - *xr2;
      *x = tmp1;
      return;
    }
    else if (logm==0) return; /* length m =1 */
  }

  /* compute a few constants */
  m=1<<logm; m2=m/2;m4=m2/2;m8=m4/2;

  /* build tables of butterfly coefficients if necessary */
  if ((logm >= 4) && (tab[logm-4]==NULL)) {
    /* allocate memory for tables */
    nel=m4-2;
    if ((tab[logm-4] = (float *) calloc(3 * nel, sizeof(float))) == NULL ) {
      fprintf(stderr,"cosine table memory error\n");
      error_exit();
    }
    /* initialize pointers */
    cn = tab[logm-4]; spcn = cn+nel; smcn =spcn + nel;

    /* compute tables */
    for (n = 1; n<m4;n++) {
      if (n==m8) continue;
      ang = n * TWOPI/m;
      c=cos(ang); s=sin(ang);
      *cn++ = c; *spcn++ = -(s+c); *smcn++ = s-c;
    }
  }

  /* step 1 */
  xr1 = x; xr2 = xr1+m2;
  for (n=0; n< m2; n++) {
    tmp1 = *xr1 + *xr2;
    *xr2 = *xr1- *xr2;
    *xr1 = tmp1;
    xr1++; xr2++;
  }

  /* step 2 */
  xr1 = x  + m2 + m4;
  for (n=0;n<m4;n++) {
    *xr1 = - *xr1;
    xr1++;
  }

  /* steps 3 and 4 */
  xr1 = x + m2; xi1 = xr1 + m4;
  if (logm >= 4) {
    nel = m4-2;
    cn = tab[logm-4]; spcn = cn + nel; smcn = spcn + nel;
  }
  xr1++; xi1++;
  for (n = 1 ; n<m4; n++) {
    if (n == m8) {
      tmp1 = SQHALF * ( *xr1 + *xi1);
      *xi1 = SQHALF * ( *xi1 - *xr1);
      *xr1 = tmp1;
    } else {
      tmp2 = *cn++ * (*xr1 + *xi1);
      tmp1 = *spcn++ * *xr1 + tmp2;
      *xr1 = *smcn++ * *xi1 + tmp2;
      *xi1 = tmp1;
    }
    xr1++; xi1++;
  }

  /* call rsrec again with half DFT length */
  rsrec(x, logm-1);
  
  /* call complex DFT routine with quarter DFT length */
  m =1 <<logm; m2 =m/2; m4 = 3 * (m/4);
  srrec(x+m2, x+m4,logm-2);

  /* step 5. sign change and data reorder */
  m = 1 <<logm; m2 = m/2; m4=m2/2; m8=m4/2;
  xr1 = x + m2 + m4;
  xr2 = x + m -1;
  for (n=0;n<m8;n++) {
    tmp1 = *xr1;
    *xr1++ = - *xr2;
    *xr2-- = - tmp1;
  }
  xr1 = x + m2 + 1;
  xr2 = x + m - 2;
  for (n=0; n<m8;n++) {
    tmp1 = *xr1;
    *xr1++ = - *xr2;
    *xr2-- = tmp1;
    xr1++;
    xr2--;
  }
  if (logm == 2) x[3]=-x[3];
}

/* direct transform for real inputs */
void rsfft(float *x, int logm) {
  /* call recursive routine */
  rsrec(x,logm);

  /* output array unshuffling using bit reversed indices */
  if (logm >1 ) {
    BR_permute (x,logm);
  }
}
/*********************************************/


