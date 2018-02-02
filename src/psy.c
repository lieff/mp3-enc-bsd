/**********************************************************************
 * ISO MPEG Audio Subgroup Software Simulation Group (1996)
 * ISO 13818-3 MPEG-2 Audio Encoder - Lower Sampling Frequency Extension
 *
 * $Id: psy.c,v 1.1 1996/02/14 04:04:23 rowlands Exp $
 *
 * $Log: psy.c,v $
 * Revision 1.1  1996/02/14 04:04:23  rowlands
 * Initial revision
 *
 * Received from Mike Coleman
 **********************************************************************/
/**********************************************************************
 *   date   programmers         comment                               *
 * 2/25/91  Davis Pan           start of version 1.0 records          *
 * 5/10/91  W. Joseph Carter    Ported to Macintosh and Unix.         *
 * 7/10/91  Earle Jennings      Ported to MsDos.                      *
 *                              replace of floats with FLOAT          *
 * 2/11/92  W. Joseph Carter    Fixed mem_alloc() arg for "absthr".   *
 * 7/24/92  M. Iwadare          HANN window coefficients modified.    *
 * 7/27/92  Masahiro Iwadare    Bug fix, FFT modification for Layer 3 *
 * 7/27/92  Masahiro Iwadare    Bug fix, "new", "old", and "oldest"   *
 *                              updates                               *
 * 8/07/92  Mike Coleman        Bug fix, read_absthr()                *
 * 95/3/21  Jon Rowlands        Removed extra debug statements        *
 **********************************************************************/

#include "common.h"
#include "encoder.h"

extern float absthr_0[];
extern float absthr_1[];
extern float absthr_2[];

FILE		*fpo;	/* file pointer */

void psycho_anal(buffer,savebuf,chn,lay,snr32,sfreq)
short int *buffer;
short int savebuf[1056];
int   chn, lay;
FLOAT snr32[32];
double sfreq;        /* to match prototype : float args are always double */
{
 unsigned int   i, j, k;
 FLOAT          r_prime, phi_prime;
 FLOAT          freq_mult, bval_lo, minthres, sum_energy;
 double         tb, temp1, temp2, temp3;

/* The static variables "r", "phi_sav", "new", "old" and "oldest" have    */
/* to be remembered for the unpredictability measure.  For "r" and        */
/* "phi_sav", the first index from the left is the channel select and     */
/* the second index is the "age" of the data.                             */

 static int     new = 0, old = 1, oldest = 0;
 static int     init = 0, flush, sync_flush, syncsize, sfreq_idx;

/* The following static variables are constants.                           */

 static double  nmt = 5.5;

 static FLOAT   crit_band[27] = {0,  100,  200, 300, 400, 510, 630,  770,
                               920, 1080, 1270,1480,1720,2000,2320, 2700,
                              3150, 3700, 4400,5300,6400,7700,9500,12000,
                             15500,25000,30000};

 static FLOAT   bmax[27] = {20.0, 20.0, 20.0, 20.0, 20.0, 17.0, 15.0,
                            10.0,  7.0,  4.4,  4.5,  4.5,  4.5,  4.5,
                             4.5,  4.5,  4.5,  4.5,  4.5,  4.5,  4.5,
                             4.5,  4.5,  4.5,  3.5,  3.5,  3.5};

/* The following pointer variables point to large areas of memory         */
/* dynamically allocated by the mem_alloc() function.  Dynamic memory     */
/* allocation is used in order to avoid stack frame or data area          */
/* overflow errors that otherwise would have occurred at compile time     */
/* on the Macintosh computer.                                             */

 FLOAT          *grouped_c, *grouped_e, *nb, *cb, *ecb, *bc;
 FLOAT          *wsamp_r, *wsamp_i, *phi, *energy;
 FLOAT          *c, *fthr;
 F32            *snrtmp;

 static int     *numlines;
 static int     *partition;
 static FLOAT   *cbval, *rnorm;
 static FLOAT   *window;
 static FLOAT   *absthr;
 static double  *tmn;
 static FCB     *s;
 static FHBLK   *lthr;
 static F2HBLK  *r, *phi_sav;

/* These dynamic memory allocations simulate "automatic" variables        */
/* placed on the stack.  For each mem_alloc() call here, there must be    */
/* a corresponding mem_free() call at the end of this function.           */

 grouped_c = (FLOAT *) mem_alloc(sizeof(FCB), "grouped_c");
 grouped_e = (FLOAT *) mem_alloc(sizeof(FCB), "grouped_e");
 nb = (FLOAT *) mem_alloc(sizeof(FCB), "nb");
 cb = (FLOAT *) mem_alloc(sizeof(FCB), "cb");
 ecb = (FLOAT *) mem_alloc(sizeof(FCB), "ecb");
 bc = (FLOAT *) mem_alloc(sizeof(FCB), "bc");
 wsamp_r = (FLOAT *) mem_alloc(sizeof(FBLK), "wsamp_r");
 wsamp_i = (FLOAT *) mem_alloc(sizeof(FBLK), "wsamp_i");
 phi = (FLOAT *) mem_alloc(sizeof(FBLK), "phi");
 energy = (FLOAT *) mem_alloc(sizeof(FBLK), "energy");
 c = (FLOAT *) mem_alloc(sizeof(FHBLK), "c");
 fthr = (FLOAT *) mem_alloc(sizeof(FHBLK), "fthr");
 snrtmp = (F32 *) mem_alloc(sizeof(F2_32), "snrtmp");

 if(init==0){

/* These dynamic memory allocations simulate "static" variables placed    */
/* in the data space.  Each mem_alloc() call here occurs only once at     */
/* initialization time.  The mem_free() function must not be called.      */

     numlines = (int *) mem_alloc(sizeof(ICB), "numlines");
     partition = (int *) mem_alloc(sizeof(IHBLK), "partition");
     fpo = fopen("out.dat", "wb");
	if(fpo==NULL) {
		puts("\t The attempt to open the output file failed.\n");
		exit(-1);}
     cbval = (FLOAT *) mem_alloc(sizeof(FCB), "cbval");
     rnorm = (FLOAT *) mem_alloc(sizeof(FCB), "rnorm");
     window = (FLOAT *) mem_alloc(sizeof(FBLK), "window");
     tmn = (double *) mem_alloc(sizeof(DCB), "tmn");
     s = (FCB *) mem_alloc(sizeof(FCBCB), "s");
     lthr = (FHBLK *) mem_alloc(sizeof(F2HBLK), "lthr");
     r = (F2HBLK *) mem_alloc(sizeof(F22HBLK), "r");
     phi_sav = (F2HBLK *) mem_alloc(sizeof(F22HBLK), "phi_sav");

     i = sfreq + 0.5;
     switch(i){
        case 32000: sfreq_idx = 0; break;
        case 44100: sfreq_idx = 1; break;
        case 48000: sfreq_idx = 2; break;
        default:    printf("error, invalid sampling frequency: %d Hz\n",i);
        exit(-1);
     }
     
     switch (sfreq_idx) {
     case 0: absthr = absthr_0; break;
     case 1: absthr = absthr_1; break;
     case 2: absthr = absthr_2; break;
     }
     
     if(lay==1){
        flush = 384;
        syncsize = 1024;
        sync_flush = 576;
     }
     else {
        flush = 384*3.0/2.0;
        syncsize = 1056;
        sync_flush = syncsize - flush;
     }
/* calculate HANN window coefficients */
/*   for(i=0;i<BLKSIZE;i++)window[i]=0.5*(1-cos(2.0*PI*i/(BLKSIZE-1.0))); */
     for(i=0;i<BLKSIZE;i++)window[i]=0.5*(1-cos(2.0*PI*(i-0.5)/BLKSIZE));
/* reset states used in unpredictability measure */
     for(i=0;i<HBLKSIZE;i++){
        r[0][0][i]=r[1][0][i]=r[0][1][i]=r[1][1][i]=0;
        phi_sav[0][0][i]=phi_sav[1][0][i]=0;
        phi_sav[0][1][i]=phi_sav[1][1][i]=0;
        lthr[0][i] = 60802371420160.0;
        lthr[1][i] = 60802371420160.0;
     }
/*****************************************************************************
 * Initialization: Compute the following constants for use later             *
 *    partition[HBLKSIZE] = the partition number associated with each        *
 *                          frequency line                                   *
 *    cbval[CBANDS]       = the center (average) bark value of each          *
 *                          partition                                        *
 *    numlines[CBANDS]    = the number of frequency lines in each partition  *
 *    tmn[CBANDS]         = tone masking noise                               *
 *****************************************************************************/
/* compute fft frequency multiplicand */
     freq_mult = sfreq/BLKSIZE;
 
/* calculate fft frequency, then bval of each line (use fthr[] as tmp storage)*/
     for(i=0;i<HBLKSIZE;i++){
        temp1 = i*freq_mult;
        j = 1;
        while(temp1>crit_band[j])j++;
        fthr[i]=j-1+(temp1-crit_band[j-1])/(crit_band[j]-crit_band[j-1]);
     }
     partition[0] = 0;
/* temp2 is the counter of the number of frequency lines in each partition */
     temp2 = 1;
     cbval[0]=fthr[0];
     bval_lo=fthr[0];
     for(i=1;i<HBLKSIZE;i++){
        if((fthr[i]-bval_lo)>0.33){
           partition[i]=partition[i-1]+1;
           cbval[partition[i-1]] = cbval[partition[i-1]]/temp2;
           cbval[partition[i]] = fthr[i];
           bval_lo = fthr[i];
           numlines[partition[i-1]] = temp2;
           temp2 = 1;
        }
        else {
           partition[i]=partition[i-1];
           cbval[partition[i]] += fthr[i];
           temp2++;
        }
     }
     numlines[partition[i-1]] = temp2;
     cbval[partition[i-1]] = cbval[partition[i-1]]/temp2;
 
/************************************************************************
 * Now compute the spreading function, s[j][i], the value of the spread-*
 * ing function, centered at band j, for band i, store for later use    *
 ************************************************************************/
     for(j=0;j<CBANDS;j++){
        for(i=0;i<CBANDS;i++){
           temp1 = (cbval[i] - cbval[j])*1.05;
           if(temp1>=0.5 && temp1<=2.5){
              temp2 = temp1 - 0.5;
              temp2 = 8.0 * (temp2*temp2 - 2.0 * temp2);
           }
           else temp2 = 0;
           temp1 += 0.474;
           temp3 = 15.811389+7.5*temp1-17.5*sqrt((double) (1.0+temp1*temp1));
           if(temp3 <= -100) s[i][j] = 0;
           else {
              temp3 = (temp2 + temp3)*LN_TO_LOG10;
              s[i][j] = exp(temp3);
           }
        }
     }

  /* Calculate Tone Masking Noise values */
     for(j=0;j<CBANDS;j++){
        temp1 = 15.5 + cbval[j];
        tmn[j] = (temp1>24.5) ? temp1 : 24.5;
  /* Calculate normalization factors for the net spreading functions */
        rnorm[j] = 0;
        for(i=0;i<CBANDS;i++){
           rnorm[j] += s[j][i];
        }
     }
     init++;
 }
 
/************************* End of Initialization *****************************/
 switch(lay) {
  case 1:
  case 2:
     for(i=0; i<lay; i++){
/*****************************************************************************
 * Net offset is 480 samples (1056-576) for layer 2; this is because one must*
 * stagger input data by 256 samples to synchronize psychoacoustic model with*
 * filter bank outputs, then stagger so that center of 1024 FFT window lines *
 * up with center of 576 "new" audio samples.                                *
 *                                                                           *
 * For layer 1, the input data still needs to be staggered by 256 samples,   *
 * then it must be staggered again so that the 384 "new" samples are centered*
 * in the 1024 FFT window.  The net offset is then 576 and you need 448 "new"*
 * samples for each iteration to keep the 384 samples of interest centered   *
 *****************************************************************************/
        for(j=0; j<syncsize; j++){
           if(j<(sync_flush))savebuf[j] = savebuf[j+flush];
           else savebuf[j] = *buffer++;
           if(j<BLKSIZE){
/**window data with HANN window***********************************************/
              wsamp_r[j] = window[j]*((FLOAT) savebuf[j]);
              wsamp_i[j] = 0;
           }
        }
/**Compute FFT****************************************************************/
        fft(wsamp_r,wsamp_i,energy,phi,1024, 0);
/*****************************************************************************
 * calculate the unpredictability measure, given energy[f] and phi[f]        *
 *****************************************************************************/
/*only update data "age" pointers after you are done with both channels      */
/*for layer 1 computations, for the layer 2 double computations, the pointers*/
/*are reset automatically on the second pass                                 */
         if(lay==2 || (lay==1 && chn==0) ){
           if(new==0){new = 1; oldest = 1;}
           else {new = 0; oldest = 0;}
           if(old==0)old = 1; else old = 0;
        }
        for(j=0; j<HBLKSIZE; j++){
           r_prime = 2.0 * r[chn][old][j] - r[chn][oldest][j];
           phi_prime = 2.0 * phi_sav[chn][old][j] - phi_sav[chn][oldest][j];
           r[chn][new][j] = sqrt((double) energy[j]);
           phi_sav[chn][new][j] = phi[j];
temp1=r[chn][new][j] * cos((double) phi[j]) - r_prime * cos((double) phi_prime);
temp2=r[chn][new][j] * sin((double) phi[j]) - r_prime * sin((double) phi_prime);
           temp3=r[chn][new][j] + fabs((double)r_prime);
           if(temp3 != 0)c[j]=sqrt(temp1*temp1+temp2*temp2)/temp3;
           else c[j] = 0;
        }
/*****************************************************************************
 * Calculate the grouped, energy-weighted, unpredictability measure,         *
 * grouped_c[], and the grouped energy. grouped_e[]                          *
 *****************************************************************************/
        for(j=1;j<CBANDS;j++){
           grouped_e[j] = 0;
           grouped_c[j] = 0;
        }
        grouped_e[0] = energy[0];
        grouped_c[0] = energy[0]*c[0];
        for(j=1;j<HBLKSIZE;j++){
           grouped_e[partition[j]] += energy[j];
           grouped_c[partition[j]] += energy[j]*c[j];
        }

/*****************************************************************************
 * convolve the grouped energy-weighted unpredictability measure             *
 * and the grouped energy with the spreading function, s[j][k]               *
 *****************************************************************************/
        for(j=0;j<CBANDS;j++){
           ecb[j] = 0;
           cb[j] = 0;
           for(k=0;k<CBANDS;k++){
              if(s[j][k] != 0.0){
                 ecb[j] += s[j][k]*grouped_e[k];
                 cb[j] += s[j][k]*grouped_c[k];
              }
           }
           if(ecb[j] !=0)cb[j] = cb[j]/ecb[j];
           else cb[j] = 0;
        }

/*****************************************************************************
 * Calculate the required SNR for each of the frequency partitions           *
 *         this whole section can be accomplished by a table lookup          *
 *****************************************************************************/
        for(j=0;j<CBANDS;j++){
           if(cb[j]<.05)cb[j]=0.05;
           else if(cb[j]>.5)cb[j]=0.5;
           tb = -0.434294482*log((double) cb[j])-0.301029996;
	   cb[j]=tb;
           bc[j] = tmn[j]*tb + nmt*(1.0-tb);
           k = cbval[j] + 0.5;
           bc[j] = (bc[j] > bmax[k]) ? bc[j] : bmax[k];
           bc[j] = exp((double) -bc[j]*LN_TO_LOG10);
        }

/*****************************************************************************
 * Calculate the permissible noise energy level in each of the frequency     *
 * partitions. Include absolute threshold and pre-echo controls              *
 *         this whole section can be accomplished by a table lookup          *
 *****************************************************************************/
        for(j=0;j<CBANDS;j++)
           if(rnorm[j] && numlines[j])
              nb[j] = ecb[j]*bc[j]/(rnorm[j]*numlines[j]);
           else nb[j] = 0;
        for(j=0;j<HBLKSIZE;j++){
/*temp1 is the preliminary threshold */
           temp1=nb[partition[j]];
           temp1=(temp1>absthr[j])?temp1:absthr[j];
/*do not use pre-echo control for layer 2 because it may do bad things to the*/
/*  MUSICAM bit allocation algorithm                                         */
           if(lay==1){
              fthr[j] = (temp1 < lthr[chn][j]) ? temp1 : lthr[chn][j];
              temp2 = temp1 * 0.00316;
              fthr[j] = (temp2 > fthr[j]) ? temp2 : fthr[j];
           }
           else fthr[j] = temp1;
           lthr[chn][j] = LXMIN*temp1;
        }

/*****************************************************************************
 * Translate the 512 threshold values to the 32 filter bands of the coder    *
 *****************************************************************************/
        for(j=0;j<193;j += 16){
           minthres = 60802371420160.0;
           sum_energy = 0.0;
           for(k=0;k<17;k++){
              if(minthres>fthr[j+k])minthres = fthr[j+k];
              sum_energy += energy[j+k];
           }
           snrtmp[i][j/16] = sum_energy/(minthres * 17.0);
           snrtmp[i][j/16] = 4.342944819 * log((double)snrtmp[i][j/16]);
        }
        for(j=208;j<(HBLKSIZE-1);j += 16){
           minthres = 0.0;
           sum_energy = 0.0;
           for(k=0;k<17;k++){
              minthres += fthr[j+k];
              sum_energy += energy[j+k];
           }
           snrtmp[i][j/16] = sum_energy/minthres;
           snrtmp[i][j/16] = 4.342944819 * log((double)snrtmp[i][j/16]);
        }
/*****************************************************************************
 * End of Psychoacuostic calculation loop                                    *
 *****************************************************************************/
     }
     for(i=0; i<32; i++){
        if(lay==2)
           snr32[i]=(snrtmp[0][i]>snrtmp[1][i])?snrtmp[0][i]:snrtmp[1][i];
        else snr32[i]=snrtmp[0][i];
     }
     break;
  case 3:
     printf("layer 3 is not currently supported\n");
     break;
  default:
     printf("error, invalid MPEG/audio coding layer: %d\n",lay);
 }

/* These mem_free() calls must correspond with the mem_alloc() calls     */
/* used at the beginning of this function to simulate "automatic"        */
/* variables placed on the stack.                                        */

 mem_free((void **) &grouped_c);
 mem_free((void **) &grouped_e);
 mem_free((void **) &nb);
 mem_free((void **) &cb);
 mem_free((void **) &ecb);
 mem_free((void **) &bc);
 mem_free((void **) &wsamp_r);
 mem_free((void **) &wsamp_i);
 mem_free((void **) &phi);
 mem_free((void **) &energy);
 mem_free((void **) &c);
 mem_free((void **) &fthr);
 mem_free((void **) &snrtmp);
}
