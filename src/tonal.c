/**********************************************************************
 * ISO MPEG Audio Subgroup Software Simulation Group (1996)
 * ISO 13818-3 MPEG-2 Audio Encoder - Lower Sampling Frequency Extension
 *
 * $Id: tonal.c,v 1.1 1996/02/14 04:04:23 rowlands Exp $
 *
 * $Log: tonal.c,v $
 * Revision 1.1  1996/02/14 04:04:23  rowlands
 * Initial revision
 *
 * Received from Mike Coleman
 **********************************************************************/
/**********************************************************************
 *   date   programmers         comment                               *
 * 2/25/91  Douglas Wong        start of version 1.1 records          *
 * 3/06/91  Douglas Wong        rename: setup.h to endef.h            *
 *                              updated I_psycho_one and II_psycho_one*
 * 3/11/91  W. J. Carter        Added Douglas Wong's updates dated    *
 *                              3/9/91 for I_Psycho_One() and for     *
 *                              II_Psycho_One().                      *
 * 5/10/91  W. Joseph Carter    Ported to Macintosh and Unix.         *
 *                              Located and fixed numerous software   *
 *                              bugs and table data errors.           *
 * 6/11/91  Davis Pan           corrected several bugs                *
 *                              based on comments from H. Fuchs       *
 * 01jul91  dpwe (Aware Inc.)   Made pow() args float                 *
 *                              Removed logical bug in I_tonal_label: *
 *                              Sometimes *tone returned == STOP      *
 * 7/10/91  Earle Jennings      no change necessary in port to MsDos  *
 * 11sep91  dpwe@aware.com      Subtracted 90.3dB from II_f_f_t peaks *
 * 10/1/91  Peter W. Farrett    Updated II_Psycho_One(),I_Psycho_One()*
 *                              to include comments.                  *
 *11/29/91  Masahiro Iwadare    Bug fix regarding POWERNORM           *
 *                              fixed several other miscellaneous bugs*
 * 2/11/92  W. Joseph Carter    Ported new code to Macintosh.  Most   *
 *                              important fixes involved changing     *
 *                              16-bit ints to long or unsigned in    *
 *                              bit alloc routines for quant of 65535 *
 *                              and passing proper function args.     *
 *                              Removed "Other Joint Stereo" option   *
 *                              and made bitrate be total channel     *
 *                              bitrate, irrespective of the mode.    *
 *                              Fixed many small bugs & reorganized.  *
 * 2/12/92  Masahiro Iwadare    Fixed some potential bugs in          *
 *          Davis Pan           subsampling()                         *
 * 2/25/92  Masahiro Iwadare    Fixed some more potential bugs        *
 * 6/24/92  Tan Ah Peng         Modified window for FFT               * 
 *                              (denominator N-1 to N)                *
 *                              Updated all critical band rate &      *
 *                              absolute threshold tables and critical*
 *                              boundaries for use with Layer I & II  *  
 *                              Corrected boundary limits for tonal   *
 *                              component computation                 *
 *                              Placement of non-tonal component at   *
 *                              geometric mean of critical band       *
 *                              (previous placement method commented  *
 *                               out - can be used if desired)        *
 * 3/01/93  Mike Li             Infinite looping fix in noise_label() *
 * 3/19/93  Jens Spille         fixed integer overflow problem in     *
 *                              psychoacoutic model 1                 *
 * 3/19/93  Giorgio Dimino      modifications to better account for   *
 *                              tonal and non-tonal components        *
 * 5/28/93 Sriram Jayasimha     "London" mod. to psychoacoustic model1*
 * 8/05/93 Masahiro Iwadare     noise_label modification "option"     *
 * 1/21/94 Seymore Shlien       fixed another infinite looping problem*
 * 7/12/95 Soeren H. Nielsen    Changes for LSF, new tables           *
 **********************************************************************/

#include "common.h"
#include "encoder.h"
#define LONDON                  /* enable "LONDON" modification */
#define MAKE_SENSE              /* enable "MAKE_SENSE" modification */
#define MI_OPTION               /* enable "MI_OPTION" modification */
/**********************************************************************
*
*        This module implements the psychoacoustic model I for the
* MPEG encoder layer II. It uses simplified tonal and noise masking
* threshold analysis to generate SMR for the encoder bit allocation
* routine.
*
**********************************************************************/

int crit_band;
int FAR *cbound;
int sub_size;

void read_cbound(lay,freq)  /* this function reads in critical */
int lay, freq;              /* band boundaries                 */
{
 int i,j,k;
 FILE *fp;
 char r[16], t[80];

 strcpy(r, "2cb1");
 r[0] = (char) lay + '0';
 r[3] = (char) freq + '0';
 if( !(fp = OpenTableFile(r)) ){       /* check boundary values */
    printf("Please check %s boundary table\n",r);
    exit(1);
 }
 fgets(t,80,fp);               /* read input for critical bands */
 sscanf(t,"%d\n",&crit_band);
 cbound = (int FAR *) mem_alloc(sizeof(int) * crit_band, "cbound");
 for(i=0;i<crit_band;i++){   /* continue to read input for */
    fgets(t,80,fp);            /* critical band boundaries   */
    sscanf(t,"%d %d\n",&j, &k);
    if(i==j) cbound[j] = k;
    else {                     /* error */
       printf("Please check index %d in cbound table %s\n",i,r);
       exit(1);
    }
 }
 fclose(fp);
}        

void read_freq_band(ltg,lay,freq)  /* this function reads in   */
int lay, freq;                     /* frequency bands and bark */
g_ptr FAR *ltg;                /* values                   */
{
 int i,j, k;
 double b,c;
 FILE *fp;
 char r[16], t[80];

 strcpy(r, "2th1");
 r[0] = (char) lay + '0';
 r[3] = (char) freq + '0';
 if( !(fp = OpenTableFile(r)) ){   /* check freq. values  */
    printf("Please check frequency and cband table %s\n",r);
    exit(1);
 }
 fgets(t,80,fp);              /* read input for freq. subbands */
 sscanf(t,"%d\n",&sub_size);
 *ltg = (g_ptr FAR ) mem_alloc(sizeof(g_thres) * sub_size, "ltg");
 (*ltg)[0].line = 0;          /* initialize global masking threshold */
 (*ltg)[0].bark = 0;
 (*ltg)[0].hear = 0;
 for(i=1;i<sub_size;i++){    /* continue to read freq. subband */
    fgets(t,80,fp);          /* and assign                     */
    sscanf(t,"%d %d %lf %lf\n",&j, &k, &b, &c);
    if(i == j){
       (*ltg)[j].line = k;
       (*ltg)[j].bark = b;
       (*ltg)[j].hear = c;
    }
    else {                   /* error */
       printf("Please check index %d in freq-cb table %s\n",i,r);
       exit(1);
    }
 }
 fclose(fp);
}

void make_map(power, ltg)       /* this function calculates the */
mask FAR power[HAN_SIZE];   /* global masking threshold     */
g_thres FAR *ltg;
{
 int i,j;

 for(i=1;i<sub_size;i++) for(j=ltg[i-1].line;j<=ltg[i].line;j++)
    power[j].map = i;
}

double add_db(a,b)
double a,b;
{
 a = pow(10.0,a/10.0);
 b = pow(10.0,b/10.0);
 return 10 * log10(a+b);
}

/****************************************************************
*
*        Fast Fourier transform of the input samples.
*
****************************************************************/

void II_f_f_t(sample, power)      /* this function calculates an */
double FAR sample[FFT_SIZE];  /* FFT analysis for the freq.  */
mask FAR power[HAN_SIZE];     /* domain                      */
{
 int i,j,k,L,l=0;
 int ip, le, le1;
 double t_r, t_i, u_r, u_i;
 static int M, MM1, init = 0, N;
 double *x_r, *x_i, *energy;
 static int *rev;
 static double *w_r, *w_i;

 x_r = (double *) mem_alloc(sizeof(DFFT), "x_r");
 x_i = (double *) mem_alloc(sizeof(DFFT), "x_i");
 energy = (double *) mem_alloc(sizeof(DFFT), "energy");
 for(i=0;i<FFT_SIZE;i++) x_r[i] = x_i[i] = energy[i] = 0;
 if(!init){
    rev = (int *) mem_alloc(sizeof(IFFT), "rev");
    w_r = (double *) mem_alloc(sizeof(D10), "w_r");
    w_i = (double *) mem_alloc(sizeof(D10), "w_i");
    M = 10;
    MM1 = 9;
    N = FFT_SIZE;
    for(L=0;L<M;L++){
       le = 1 << (M-L);
       le1 = le >> 1;
       w_r[L] = cos(PI/le1);
       w_i[L] = -sin(PI/le1);
    }
    for(i=0;i<FFT_SIZE;rev[i] = l,i++) for(j=0,l=0;j<10;j++){
       k=(i>>j) & 1;
       l |= (k<<(9-j));                
    }
    init = 1;
 }
 memcpy( (char *) x_r, (char *) sample, sizeof(double) * FFT_SIZE);
 for(L=0;L<MM1;L++){
    le = 1 << (M-L);
    le1 = le >> 1;
    u_r = 1;
    u_i = 0;
    for(j=0;j<le1;j++){
       for(i=j;i<N;i+=le){
          ip = i + le1;
          t_r = x_r[i] + x_r[ip];
          t_i = x_i[i] + x_i[ip];
          x_r[ip] = x_r[i] - x_r[ip];
          x_i[ip] = x_i[i] - x_i[ip];
          x_r[i] = t_r;
          x_i[i] = t_i;
          t_r = x_r[ip];
          x_r[ip] = x_r[ip] * u_r - x_i[ip] * u_i;
          x_i[ip] = x_i[ip] * u_r + t_r * u_i;
       }
       t_r = u_r;
       u_r = u_r * w_r[L] - u_i * w_i[L];
       u_i = u_i * w_r[L] + t_r * w_i[L];
    }
 }
 for(i=0;i<N;i+=2){
    ip = i + 1;
    t_r = x_r[i] + x_r[ip];
    t_i = x_i[i] + x_i[ip];
    x_r[ip] = x_r[i] - x_r[ip];
    x_i[ip] = x_i[i] - x_i[ip];
    x_r[i] = t_r;
    x_i[i] = t_i;
    energy[i] = x_r[i] * x_r[i] + x_i[i] * x_i[i];
 }
 for(i=0;i<FFT_SIZE;i++) if(i<rev[i]){
    t_r = energy[i];
    energy[i] = energy[rev[i]];
    energy[rev[i]] = t_r;
 }
 for(i=0;i<HAN_SIZE;i++){    /* calculate power density spectrum */
    if (energy[i] < 1E-20) energy[i] = 1E-20;
    power[i].x = 10 * log10(energy[i]) + POWERNORM;
    power[i].next = STOP;
    power[i].type = FALSE;
 }
 mem_free((void **) &x_r);
 mem_free((void **) &x_i);
 mem_free((void **) &energy);
}

/****************************************************************
*
*         Window the incoming audio signal.
*
****************************************************************/

void II_hann_win(sample)          /* this function calculates a  */
double FAR sample[FFT_SIZE];  /* Hann window for PCM (input) */
{                                 /* samples for a 1024-pt. FFT  */
 register int i;
 register double sqrt_8_over_3;
 static int init = 0;
 static double FAR *window;

 if(!init){  /* calculate window function for the Fourier transform */
    window = (double FAR *) mem_alloc(sizeof(DFFT), "window");
    sqrt_8_over_3 = pow(8.0/3.0, 0.5);
    for(i=0;i<FFT_SIZE;i++){
       /* Hann window formula */
       window[i]=sqrt_8_over_3*0.5*(1-cos(2.0*PI*i/(FFT_SIZE)))/FFT_SIZE;
    }
    init = 1;
 }
 for(i=0;i<FFT_SIZE;i++) sample[i] *= window[i];
}

/*******************************************************************
*
*        This function finds the maximum spectral component in each
* subband and return them to the encoder for time-domain threshold
* determination.
*
*******************************************************************/
#ifndef LONDON
void II_pick_max(power, spike)
double FAR spike[SBLIMIT];
mask FAR power[HAN_SIZE];
{
 double max;
 int i,j;

 for(i=0;i<HAN_SIZE;spike[i>>4] = max, i+=16)      /* calculate the      */
 for(j=0, max = DBMIN;j<16;j++)                    /* maximum spectral   */
    max = (max>power[i+j].x) ? max : power[i+j].x; /* component in each  */
}                                                  /* subband from bound */
                                                   /* 4-16               */
#else
void II_pick_max(power, spike)
double FAR spike[SBLIMIT];
mask FAR power[HAN_SIZE];
{
 double sum;
 int i,j;

 for(i=0;i<HAN_SIZE;spike[i>>4] = 10.0*log10(sum), i+=16)
                                                   /* calculate the      */
 for(j=0, sum = pow(10.0,0.1*DBMIN);j<16;j++)      /* sum of spectral   */
   sum += pow(10.0,0.1*power[i+j].x);              /* component in each  */
}                                                  /* subband from bound */
                                                   /* 4-16               */
#endif

/****************************************************************
*
*        This function labels the tonal component in the power
* spectrum.
*
****************************************************************/

void II_tonal_label(power, tone)  /* this function extracts (tonal) */
mask FAR power[HAN_SIZE];     /* sinusoidals from the spectrum  */
int *tone;
{
 int i,j, last = LAST, first, run, last_but_one = LAST; /* dpwe */
 double max;

 *tone = LAST;
 for(i=2;i<HAN_SIZE-12;i++){
    if(power[i].x>power[i-1].x && power[i].x>=power[i+1].x){
       power[i].type = TONE;
       power[i].next = LAST;
       if(last != LAST) power[last].next = i;
       else first = *tone = i;
       last = i;
    }
 }
 last = LAST;
 first = *tone;
 *tone = LAST;
 while(first != LAST){               /* the conditions for the tonal          */
    if(first<3 || first>500) run = 0;/* otherwise k+/-j will be out of bounds */
    else if(first<63) run = 2;       /* components in layer II, which         */
    else if(first<127) run = 3;      /* are the boundaries for calc.          */
    else if(first<255) run = 6;      /* the tonal components                  */
    else run = 12;
    max = power[first].x - 7;        /* after calculation of tonal   */
    for(j=2;j<=run;j++)              /* components, set to local max */
       if(max < power[first-j].x || max < power[first+j].x){
          power[first].type = FALSE;
          break;
       }
    if(power[first].type == TONE){   /* extract tonal components */
       int help=first;
       if(*tone==LAST) *tone = first;
       while((power[help].next!=LAST)&&(power[help].next-first)<=run)
          help=power[help].next;
       help=power[help].next;
       power[first].next=help;
       if((first-last)<=run){
          if(last_but_one != LAST) power[last_but_one].next=first;
       }
       if(first>1 && first<500){     /* calculate the sum of the */
          double tmp;                /* powers of the components */
          tmp = add_db(power[first-1].x, power[first+1].x);
          power[first].x = add_db(power[first].x, tmp);
       }
       for(j=1;j<=run;j++){
          power[first-j].x = power[first+j].x = DBMIN;
          power[first-j].next = power[first+j].next = STOP;
          power[first-j].type = power[first+j].type = FALSE;
       }
       last_but_one=last;
       last = first;
       first = power[first].next;
    }
    else {
       int ll;
       if(last == LAST); /* *tone = power[first].next; dpwe */
       else power[last].next = power[first].next;
       ll = first;
       first = power[first].next;
       power[ll].next = STOP;
    }
 }
}

/****************************************************************
*
*        This function groups all the remaining non-tonal
* spectral lines into critical band where they are replaced by
* one single line.
*
****************************************************************/
        
void noise_label(power, noise, ltg)
g_thres FAR *ltg;
mask FAR *power;
int *noise;
{
 int i,j, centre, last = LAST;
 double index, weight, sum;
                              /* calculate the remaining spectral */
 for(i=0;i<crit_band-1;i++){  /* lines for non-tonal components   */
     for(j=cbound[i],weight = 0.0,sum = DBMIN;j<cbound[i+1];j++){
        if(power[j].type != TONE){
           if(power[j].x != DBMIN){
              sum = add_db(power[j].x,sum);
/* the line below and others under the "MAKE_SENSE" condition are an alternate
   interpretation of "geometric mean". This approach may make more sense but
   it has not been tested with hardware. */
#ifdef MAKE_SENSE
/* weight += pow(10.0, power[j].x/10.0) * (ltg[power[j].map].bark-i);
   bad code [SS] 21-1-93
 */
    weight += pow(10.0,power[j].x/10.0) * (double) (j-cbound[i]) /
     (double) (cbound[i+1]-cbound[i]);  /* correction */
#endif
              power[j].x = DBMIN;
           }
        }   /*  check to see if the spectral line is low dB, and if  */
     }      /* so replace the center of the critical band, which is */
            /* the center freq. of the noise component              */

#ifdef MAKE_SENSE
     if(sum <= DBMIN)  centre = (cbound[i+1]+cbound[i]) /2;
     else {
        index = weight/pow(10.0,sum/10.0);
        centre = cbound[i] + (int) (index * (double) (cbound[i+1]-cbound[i]) );
     } 
#else
     index = (double)( ((double)cbound[i]) * ((double)(cbound[i+1]-1)) );
     centre = (int)(pow(index,0.5)+0.5);
#endif

    /* locate next non-tonal component until finished; */
    /* add to list of non-tonal components             */
#ifdef MI_OPTION
     /* Masahiro Iwadare's fix for infinite looping problem? */
     if(power[centre].type == TONE) 
       if (power[centre+1].type == TONE) centre++; else centre--;
#else
     /* Mike Li's fix for infinite looping problem */
     if(power[centre].type == FALSE) centre++;

     if(power[centre].type == NOISE){
       if(power[centre].x >= ltg[power[i].map].hear){
         if(sum >= ltg[power[i].map].hear) sum = add_db(power[j].x,sum);
         else
         sum = power[centre].x;
       }
     }
#endif
     if(last == LAST) *noise = centre;
     else {
        power[centre].next = LAST;
        power[last].next = centre;
     }
     power[centre].x = sum;
     power[centre].type = NOISE;        
     last = centre;
 }        
}

/****************************************************************
*
*        This function reduces the number of noise and tonal
* component for further threshold analysis.
*
****************************************************************/

void subsampling(power, ltg, tone, noise)
mask FAR power[HAN_SIZE];
g_thres FAR *ltg;
int *tone, *noise;
{
 int i, old;

 i = *tone; old = STOP;    /* calculate tonal components for */
 while(i!=LAST){           /* reduction of spectral lines    */
    if(power[i].x < ltg[power[i].map].hear){
       power[i].type = FALSE;
       power[i].x = DBMIN;
       if(old == STOP) *tone = power[i].next;
       else power[old].next = power[i].next;
    }
    else old = i;
    i = power[i].next;
 }
 i = *noise; old = STOP;    /* calculate non-tonal components for */
 while(i!=LAST){            /* reduction of spectral lines        */
    if(power[i].x < ltg[power[i].map].hear){
       power[i].type = FALSE;
       power[i].x = DBMIN;
       if(old == STOP) *noise = power[i].next;
       else power[old].next = power[i].next;
    }
    else old = i;
    i = power[i].next;
 }
 i = *tone; old = STOP;
 while(i != LAST){                              /* if more than one */
    if(power[i].next == LAST)break;             /* tonal component  */
    if(ltg[power[power[i].next].map].bark -     /* is less than .5  */
       ltg[power[i].map].bark < 0.5) {          /* bark, take the   */
       if(power[power[i].next].x > power[i].x ){/* maximum          */
          if(old == STOP) *tone = power[i].next;
          else power[old].next = power[i].next;
          power[i].type = FALSE;
          power[i].x = DBMIN;
          i = power[i].next;
       }
       else {
          power[power[i].next].type = FALSE;
          power[power[i].next].x = DBMIN;
          power[i].next = power[power[i].next].next;
          old = i;
       }
    }
    else {
      old = i;
      i = power[i].next;
    }
 }
}

/****************************************************************
*
*        This function calculates the individual threshold and
* sum with the quiet threshold to find the global threshold.
*
****************************************************************/

void threshold(power, ltg, tone, noise, bit_rate)
mask FAR power[HAN_SIZE];
g_thres FAR *ltg;
int *tone, *noise, bit_rate;
{
 int k, t;
 double dz, tmps, vf;

 for(k=1;k<sub_size;k++){
    ltg[k].x = DBMIN;
    t = *tone;          /* calculate individual masking threshold for */
    while(t != LAST){   /* components in order to find the global     */
       if(ltg[k].bark-ltg[power[t].map].bark >= -3.0 && /*threshold (LTG)*/
          ltg[k].bark-ltg[power[t].map].bark <8.0){
          dz = ltg[k].bark-ltg[power[t].map].bark; /* distance of bark value*/
          tmps = -1.525-0.275*ltg[power[t].map].bark - 4.5 + power[t].x;
             /* masking function for lower & upper slopes */
          if(-3<=dz && dz<-1) vf = 17*(dz+1)-(0.4*power[t].x +6);
          else if(-1<=dz && dz<0) vf = (0.4 *power[t].x + 6) * dz;
          else if(0<=dz && dz<1) vf = (-17*dz);
          else if(1<=dz && dz<8) vf = -(dz-1) * (17-0.15 *power[t].x) - 17;
          tmps += vf;        
          ltg[k].x = add_db(ltg[k].x, tmps);
       }
       t = power[t].next;
    }

    t = *noise;        /* calculate individual masking threshold  */
    while(t != LAST){  /* for non-tonal components to find LTG    */
       if(ltg[k].bark-ltg[power[t].map].bark >= -3.0 &&
          ltg[k].bark-ltg[power[t].map].bark <8.0){
          dz = ltg[k].bark-ltg[power[t].map].bark; /* distance of bark value */
          tmps = -1.525-0.175*ltg[power[t].map].bark -0.5 + power[t].x;
             /* masking function for lower & upper slopes */
          if(-3<=dz && dz<-1) vf = 17*(dz+1)-(0.4*power[t].x +6);
          else if(-1<=dz && dz<0) vf = (0.4 *power[t].x + 6) * dz;
          else if(0<=dz && dz<1) vf = (-17*dz);
          else if(1<=dz && dz<8) vf = -(dz-1) * (17-0.15 *power[t].x) - 17;
          tmps += vf;
          ltg[k].x = add_db(ltg[k].x, tmps);
       }
       t = power[t].next;
    }
    if(bit_rate<96)ltg[k].x = add_db(ltg[k].hear, ltg[k].x);
    else ltg[k].x = add_db(ltg[k].hear-12.0, ltg[k].x);
 }
}

/****************************************************************
*
*        This function finds the minimum masking threshold and
* return the value to the encoder.
*
****************************************************************/

void II_minimum_mask(ltg,ltmin,sblimit)
g_thres FAR *ltg;
double FAR ltmin[SBLIMIT];
int sblimit;
{
 double min;
 int i,j;

 j=1;
 for(i=0;i<sblimit;i++)
    if(j>=sub_size-1)                   /* check subband limit, and       */
       ltmin[i] = ltg[sub_size-1].hear; /* calculate the minimum masking  */
    else {                              /* level of LTMIN for each subband*/
       min = ltg[j].x;
       while(ltg[j].line>>4 == i && j < sub_size){
       if(min>ltg[j].x)  min = ltg[j].x;
       j++;
    }
    ltmin[i] = min;
 }
}

/*****************************************************************
*
*        This procedure is called in musicin to pick out the
* smaller of the scalefactor or threshold.
*
*****************************************************************/

void II_smr(ltmin, spike, scale, sblimit)
double FAR spike[SBLIMIT], scale[SBLIMIT], ltmin[SBLIMIT];
int sblimit;
{
 int i;
 double max;
                
 for(i=0;i<sblimit;i++){                     /* determine the signal   */
    max = 20 * log10(scale[i] * 32768) - 10; /* level for each subband */
    if(spike[i]>max) max = spike[i];         /* for the maximum scale  */
    max -= ltmin[i];                         /* factors                */
    ltmin[i] = max;
 }
}
        
/****************************************************************
*
*        This procedure calls all the necessary functions to
* complete the psychoacoustic analysis.
*
****************************************************************/

void II_Psycho_One(buffer, scale, ltmin, fr_ps)
short FAR buffer[2][1152];
double FAR scale[2][SBLIMIT], ltmin[2][SBLIMIT];
frame_params *fr_ps;
{
 layer *info = fr_ps->header;
 int   stereo = fr_ps->stereo;
 int   sblimit = fr_ps->sblimit;
 int k,i, tone=0, noise=0;
 static char init = 0;
 static int off[2] = {256,256};
 double *sample;
 DSBL *spike;
 static D1408 *fft_buf;
 static mask_ptr FAR power;
 static g_ptr FAR ltg;

 sample = (double *) mem_alloc(sizeof(DFFT), "sample");
 spike = (DSBL *) mem_alloc(sizeof(D2SBL), "spike");
     /* call functions for critical boundaries, freq. */
 if(!init){  /* bands, bark values, and mapping */
    fft_buf = (D1408 *) mem_alloc((long) sizeof(D1408) * 2, "fft_buf");
    power = (mask_ptr FAR ) mem_alloc(sizeof(mask) * HAN_SIZE, "power");
    if (info->version == MPEG_AUDIO_ID) {
      read_cbound(info->lay, info->sampling_frequency);
      read_freq_band(&ltg, info->lay, info->sampling_frequency);
    } else {
      read_cbound(info->lay, info->sampling_frequency + 4);
      read_freq_band(&ltg, info->lay, info->sampling_frequency + 4);
    }
    make_map(power,ltg);
    for (i=0;i<1408;i++) fft_buf[0][i] = fft_buf[1][i] = 0;
    init = 1;
 }
 for(k=0;k<stereo;k++){  /* check pcm input for 3 blocks of 384 samples */
    for(i=0;i<1152;i++) fft_buf[k][(i+off[k])%1408]= (double)buffer[k][i]/SCALE;
    for(i=0;i<FFT_SIZE;i++) sample[i] = fft_buf[k][(i+1216+off[k])%1408];
    off[k] += 1152;
    off[k] %= 1408;
                            /* call functions for windowing PCM samples,*/
    II_hann_win(sample);    /* location of spectral components in each  */
    for(i=0;i<HAN_SIZE;i++) power[i].x = DBMIN;  /*subband with labeling*/
    II_f_f_t(sample, power);                     /*locate remaining non-*/
    II_pick_max(power, &spike[k][0]);            /*tonal sinusoidals,   */
    II_tonal_label(power, &tone);                /*reduce noise & tonal */
    noise_label(power, &noise, ltg);             /*components, find     */
    subsampling(power, ltg, &tone, &noise);      /*global & minimal     */
    threshold(power, ltg, &tone, &noise,         /*threshold, and sgnl- */
      bitrate[info->version][info->lay-1][info->bitrate_index]/stereo); /*to-mask ratio*/
    II_minimum_mask(ltg, &ltmin[k][0], sblimit);
    II_smr(&ltmin[k][0], &spike[k][0], &scale[k][0], sblimit);        
 }
 mem_free((void **) &sample);
 mem_free((void **) &spike);
}

/**********************************************************************
*
*        This module implements the psychoacoustic model I for the
* MPEG encoder layer I. It uses simplified tonal and noise masking
* threshold analysis to generate SMR for the encoder bit allocation
* routine.
*
**********************************************************************/

/****************************************************************
*
*        Fast Fourier transform of the input samples.
*
****************************************************************/

void I_f_f_t(sample, power)         /* this function calculates */
double FAR sample[FFT_SIZE/2];  /* an FFT analysis for the  */
mask FAR power[HAN_SIZE/2];     /* freq. domain             */
{
 int i,j,k,L,l=0;
 int ip, le, le1;
 double t_r, t_i, u_r, u_i;
 static int M, MM1, init = 0, N;
 double *x_r, *x_i, *energy;
 static int *rev;
 static double *w_r, *w_i;

 x_r = (double *) mem_alloc(sizeof(DFFT2), "x_r");
 x_i = (double *) mem_alloc(sizeof(DFFT2), "x_i");
 energy = (double *) mem_alloc(sizeof(DFFT2), "energy");
 for(i=0;i<FFT_SIZE/2;i++) x_r[i] = x_i[i] = energy[i] = 0;
 if(!init){
    rev = (int *) mem_alloc(sizeof(IFFT2), "rev");
    w_r = (double *) mem_alloc(sizeof(D9), "w_r");
    w_i = (double *) mem_alloc(sizeof(D9), "w_i");
    M = 9;
    MM1 = 8;
    N = FFT_SIZE/2;
    for(L=0;L<M;L++){
       le = 1 << (M-L);
       le1 = le >> 1;
       w_r[L] = cos(PI/le1);
       w_i[L] = -sin(PI/le1);
    }
    for(i=0;i<FFT_SIZE/2;rev[i] = l,i++) for(j=0,l=0;j<9;j++){
       k=(i>>j) & 1;
       l |= (k<<(8-j));                
    }
    init = 1;
 }
 memcpy( (char *) x_r, (char *) sample, sizeof(double) * FFT_SIZE/2);
 for(L=0;L<MM1;L++){
    le = 1 << (M-L);
    le1 = le >> 1;
    u_r = 1;
    u_i = 0;
    for(j=0;j<le1;j++){
       for(i=j;i<N;i+=le){
          ip = i + le1;
          t_r = x_r[i] + x_r[ip];
          t_i = x_i[i] + x_i[ip];
          x_r[ip] = x_r[i] - x_r[ip];
          x_i[ip] = x_i[i] - x_i[ip];
          x_r[i] = t_r;
          x_i[i] = t_i;
          t_r = x_r[ip];
          x_r[ip] = x_r[ip] * u_r - x_i[ip] * u_i;
          x_i[ip] = x_i[ip] * u_r + t_r * u_i;
       }
       t_r = u_r;
       u_r = u_r * w_r[L] - u_i * w_i[L];
       u_i = u_i * w_r[L] + t_r * w_i[L];
    }
 }
 for(i=0;i<N;i+=2){
    ip = i + 1;
    t_r = x_r[i] + x_r[ip];
    t_i = x_i[i] + x_i[ip];
    x_r[ip] = x_r[i] - x_r[ip];
    x_i[ip] = x_i[i] - x_i[ip];
    x_r[i] = t_r;
    x_i[i] = t_i;
    energy[i] = x_r[i] * x_r[i] + x_i[i] * x_i[i];
 }
 for(i=0;i<FFT_SIZE/2;i++) if(i<rev[i]){
    t_r = energy[i];
    energy[i] = energy[rev[i]];
    energy[rev[i]] = t_r;
 }
 for(i=0;i<HAN_SIZE/2;i++){                     /* calculate power  */
    if(energy[i] < 1E-20) energy[i] = 1E-20;    /* density spectrum */
       power[i].x = 10 * log10(energy[i]) + POWERNORM;
       power[i].next = STOP;
       power[i].type = FALSE;
 }
 mem_free((void **) &x_r);
 mem_free((void **) &x_i);
 mem_free((void **) &energy);
}

/****************************************************************
*
*         Window the incoming audio signal.
*
****************************************************************/

void I_hann_win(sample)             /* this function calculates a  */
double FAR sample[FFT_SIZE/2];  /* Hann window for PCM (input) */
{                                   /* samples for a 512-pt. FFT   */
 register int i;
 register double sqrt_8_over_3;
 static int init = 0;
 static double FAR *window;

 if(!init){  /* calculate window function for the Fourier transform */
    window = (double FAR *) mem_alloc(sizeof(DFFT2), "window");
    sqrt_8_over_3 = pow(8.0/3.0, 0.5);
    for(i=0;i<FFT_SIZE/2;i++){
      /* Hann window formula */
      window[i]=sqrt_8_over_3*0.5*(1-cos(2.0*PI*i/(FFT_SIZE/2)))/(FFT_SIZE/2);
    }
    init = 1;
 }
 for(i=0;i<FFT_SIZE/2;i++) sample[i] *= window[i];
}

/*******************************************************************
*
*        This function finds the maximum spectral component in each
* subband and return them to the encoder for time-domain threshold
* determination.
*
*******************************************************************/
#ifndef LONDON
void I_pick_max(power, spike)
double FAR spike[SBLIMIT];
mask FAR power[HAN_SIZE/2];
{
 double max;
 int i,j;

 /* calculate the spectral component in each subband */
 for(i=0;i<HAN_SIZE/2;spike[i>>3] = max, i+=8)
    for(j=0, max = DBMIN;j<8;j++) max = (max>power[i+j].x) ? max : power[i+j].x;
}
#else
void I_pick_max(power, spike)
double FAR spike[SBLIMIT];
mask FAR power[HAN_SIZE];
{
 double sum;
 int i,j;

 for(i=0;i<HAN_SIZE/2;spike[i>>3] = 10.0*log10(sum), i+=8)
                                                   /* calculate the      */
 for(j=0, sum = pow(10.0,0.1*DBMIN);j<8;j++)       /* sum of spectral   */
   sum += pow(10.0,0.1*power[i+j].x);              /* component in each  */
}                                                  /* subband from bound */
#endif
/****************************************************************
*
*        This function labels the tonal component in the power
* spectrum.
*
****************************************************************/

void I_tonal_label(power, tone)     /* this function extracts   */
mask FAR power[HAN_SIZE/2];     /* (tonal) sinusoidals from */
int *tone;                          /* the spectrum             */
{
 int i,j, last = LAST, first, run;
 double max;
 int last_but_one= LAST;

 *tone = LAST;
 for(i=2;i<HAN_SIZE/2-6;i++){
    if(power[i].x>power[i-1].x && power[i].x>=power[i+1].x){
       power[i].type = TONE;
       power[i].next = LAST;
       if(last != LAST) power[last].next = i;
       else first = *tone = i;
       last = i;
    }
 }
 last = LAST;
 first = *tone;
 *tone = LAST;
 while(first != LAST){                /* conditions for the tonal     */
    if(first<3 || first>250) run = 0; /* otherwise k+/-j will be out of bounds*/
    else if(first<63) run = 2;        /* components in layer I, which */
    else if(first<127) run = 3;       /* are the boundaries for calc.   */
    else run = 6;                     /* the tonal components          */
    max = power[first].x - 7;
    for(j=2;j<=run;j++)  /* after calc. of tonal components, set to loc.*/
       if(max < power[first-j].x || max < power[first+j].x){   /* max   */
          power[first].type = FALSE;
          break;
       }
    if(power[first].type == TONE){    /* extract tonal components */
       int help=first;
       if(*tone == LAST) *tone = first;
       while((power[help].next!=LAST)&&(power[help].next-first)<=run)
          help=power[help].next;
       help=power[help].next;
       power[first].next=help;
       if((first-last)<=run){
          if(last_but_one != LAST) power[last_but_one].next=first;
       }
       if(first>1 && first<255){     /* calculate the sum of the */
          double tmp;                /* powers of the components */
          tmp = add_db(power[first-1].x, power[first+1].x);
          power[first].x = add_db(power[first].x, tmp);
       }
       for(j=1;j<=run;j++){
          power[first-j].x = power[first+j].x = DBMIN;
          power[first-j].next = power[first+j].next = STOP; /*dpwe: 2nd was .x*/
          power[first-j].type = power[first+j].type = FALSE;
       }
       last_but_one=last;
       last = first;
       first = power[first].next;
    }
    else {
       int ll;
       if(last == LAST) ; /* *tone = power[first].next; dpwe */
       else power[last].next = power[first].next;
       ll = first;
       first = power[first].next;
       power[ll].next = STOP;
    }
 }
}                        
                                
/****************************************************************
*
*        This function finds the minimum masking threshold and
* return the value to the encoder.
*
****************************************************************/

void I_minimum_mask(ltg,ltmin)
g_thres FAR *ltg;
double FAR ltmin[SBLIMIT];
{
 double min;
 int i,j;

 j=1;
 for(i=0;i<SBLIMIT;i++)
    if(j>=sub_size-1)                   /* check subband limit, and       */
       ltmin[i] = ltg[sub_size-1].hear; /* calculate the minimum masking  */
    else {                              /* level of LTMIN for each subband*/
       min = ltg[j].x;
       while(ltg[j].line>>3 == i && j < sub_size){
          if (min>ltg[j].x)  min = ltg[j].x;
          j++;
       }
       ltmin[i] = min;
    }
}

/*****************************************************************
*
*        This procedure is called in musicin to pick out the
* smaller of the scalefactor or threshold.
*
*****************************************************************/

void I_smr(ltmin, spike, scale)
double FAR spike[SBLIMIT], scale[SBLIMIT], ltmin[SBLIMIT];
{
 int i;
 double max;
                
 for(i=0;i<SBLIMIT;i++){                      /* determine the signal   */
    max = 20 * log10(scale[i] * 32768) - 10;  /* level for each subband */
    if(spike[i]>max) max = spike[i];          /* for the scalefactor    */
    max -= ltmin[i];
    ltmin[i] = max;
 }
}
        
/****************************************************************
*
*        This procedure calls all the necessary functions to
* complete the psychoacoustic analysis.
*
****************************************************************/

void I_Psycho_One(buffer, scale, ltmin, fr_ps)
short FAR buffer[2][1152];
double FAR scale[2][SBLIMIT], ltmin[2][SBLIMIT];
frame_params *fr_ps;
{
 int stereo = fr_ps->stereo;
 the_layer info = fr_ps->header;
 int k,i, tone=0, noise=0;
 static char init = 0;
 static int off[2] = {256,256};
 double *sample;
 DSBL *spike;
 static D640 *fft_buf;
 static mask_ptr FAR power;
 static g_ptr FAR ltg;

 sample = (double *) mem_alloc(sizeof(DFFT2), "sample");
 spike = (DSBL *) mem_alloc(sizeof(D2SBL), "spike");
            /* call functions for critical boundaries, freq. */
 if(!init){ /* bands, bark values, and mapping              */
    fft_buf = (D640 *) mem_alloc(sizeof(D640) * 2, "fft_buf");
    power = (mask_ptr FAR ) mem_alloc(sizeof(mask) * HAN_SIZE/2, "power");
    if (info->version == MPEG_AUDIO_ID) {
      read_cbound(info->lay, info->sampling_frequency);
      read_freq_band(&ltg, info->lay, info->sampling_frequency);
    } else {
      read_cbound(info->lay, info->sampling_frequency + 4);
      read_freq_band(&ltg, info->lay, info->sampling_frequency + 4);
    }
    make_map(power,ltg);
    for(i=0;i<640;i++) fft_buf[0][i] = fft_buf[1][i] = 0;
    init = 1;
 }
 for(k=0;k<stereo;k++){    /* check PCM input for a block of */
    for(i=0;i<384;i++)     /* 384 samples for a 512-pt. FFT  */
       fft_buf[k][(i+off[k])%640]= (double) buffer[k][i]/SCALE;
    for(i=0;i<FFT_SIZE/2;i++)
       sample[i] = fft_buf[k][(i+448+off[k])%640];
    off[k] += 384;
    off[k] %= 640;
                        /* call functions for windowing PCM samples,   */
    I_hann_win(sample); /* location of spectral components in each     */
    for(i=0;i<HAN_SIZE/2;i++) power[i].x = DBMIN;   /* subband with    */
    I_f_f_t(sample, power);              /* labeling, locate remaining */
    I_pick_max(power, &spike[k][0]);     /* non-tonal sinusoidals,     */
    I_tonal_label(power, &tone);         /* reduce noise & tonal com., */
    noise_label(power, &noise, ltg);     /* find global & minimal      */
    subsampling(power, ltg, &tone, &noise);  /* threshold, and sgnl-   */
    threshold(power, ltg, &tone, &noise,     /* to-mask ratio          */
      bitrate[info->version][info->lay-1][info->bitrate_index]/stereo);
    I_minimum_mask(ltg, &ltmin[k][0]);
    I_smr(&ltmin[k][0], &spike[k][0], &scale[k][0]);        
 }
 mem_free((void **) &sample);
 mem_free((void **) &spike);
}
