/**********************************************************************
 * ISO MPEG Audio Subgroup Software Simulation Group (1996)
 * ISO 13818-3 MPEG-2 Audio Encoder - Lower Sampling Frequency Extension
 *
 * $Id: l3psy.c,v 1.2 1997/01/19 22:28:29 rowlands Exp $
 *
 * $Log: l3psy.c,v $
 * Revision 1.2  1997/01/19 22:28:29  rowlands
 * Layer 3 bug fixes from Seymour Shlien
 *
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
 * 3/16/92  Masahiro Iwadare	Modification for Layer III            *
 * 17/4/93  Masahiro Iwadare    Updated for IS Modification           *
 **********************************************************************/

#include "common.h"
#include "encoder.h"
#include "l3psy.h"
#include "l3side.h"
#include <assert.h>

#define maximum(x,y) ( (x>y) ? x : y )
#define minimum(x,y) ( (x<y) ? x : y )

extern float absthr_0[];
extern float absthr_1[];
extern float absthr_2[];
extern double psy_data[];
void sprdngf1(FLOAT *, double *);
void sprdngf2(double *, FLOAT *);
static double s3_l[CBANDS][CBANDS];

void L3para_read( double sfreq, int numlines[CBANDS], int partition_l[HBLKSIZE],
		  double minval[CBANDS], double qthr_l[CBANDS], double norm_l[CBANDS],
		  double s3_l[CBANDS][CBANDS], int partition_s[HBLKSIZE_s], double qthr_s[CBANDS_s],
		  double norm_s[CBANDS_s], double SNR_s[CBANDS_s],
		  int cbw_l[SBMAX_l], int bu_l[SBMAX_l], int bo_l[SBMAX_l],
		  double w1_l[SBMAX_l], double w2_l[SBMAX_l],
		  int cbw_s[SBMAX_s], int bu_s[SBMAX_s], int bo_s[SBMAX_s],
		  double w1_s[SBMAX_s], double w2_s[SBMAX_s] );
									
void L3psycho_anal( short int *buffer, short int savebuf[1344], int chn, int lay, FLOAT snr32[32],
		    double sfreq, double ratio_d[21], double ratio_ds[12][3],
		    double *pe, gr_info *cod_info )
{
    static double ratio[2][21];
    static double ratio_s[2][12][3];
    int blocktype;
    unsigned int   b, i, j, k;
    double         r_prime, phi_prime; /* not FLOAT */
    FLOAT          freq_mult, bval_lo, min_thres, sum_energy;
    double         tb, temp1,temp2,temp3;

    /*         nint(); Layer III */
    double   thr[CBANDS];

/* The static variables "r", "phi_sav", "new", "old" and "oldest" have    */
/* to be remembered for the unpredictability measure.  For "r" and        */
/* "phi_sav", the first index from the left is the channel select and     */
/* the second index is the "age" of the data.                             */


   static FLOAT window_s[BLKSIZE_s] ;
 static int     new = 0, old = 1, oldest = 0;
 static int     init = 0, flush, sync_flush, syncsize, sfreq_idx;
 static double 	cw[HBLKSIZE], eb[CBANDS];
 static double 	ctb[CBANDS];
 static double	SNR_l[CBANDS], SNR_s[CBANDS_s];
 static int	init_L3;
 static double	minval[CBANDS],qthr_l[CBANDS],norm_l[CBANDS];
 static double	qthr_s[CBANDS_s],norm_s[CBANDS_s];
 static double	nb_1[2][CBANDS], nb_2[2][CBANDS];

/* Scale Factor Bands */
 static int	cbw_l[SBMAX_l],bu_l[SBMAX_l],bo_l[SBMAX_l] ;
 static int	cbw_s[SBMAX_s],bu_s[SBMAX_s],bo_s[SBMAX_s] ;
 static double	w1_l[SBMAX_l], w2_l[SBMAX_l];
 static double	w1_s[SBMAX_s], w2_s[SBMAX_s];
 static double	en[SBMAX_l],   thm[SBMAX_l] ;
 static int	blocktype_old[2] ;
 int	sb,sblock;
 static int	partition_l[HBLKSIZE],partition_s[HBLKSIZE_s];


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
 static FLOAT	energy_s[3][256];
 static FLOAT phi_s[3][256] ; /* 256 samples not 129 */
 FLOAT          *c, *fthr;
 F32            *snrtmp;

 static	int	*numlines ;
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

    assert( lay == 3 );
 if(init==0){

/* These dynamic memory allocations simulate "static" variables placed    */
/* in the data space.  Each mem_alloc() call here occurs only once at     */
/* initialization time.  The mem_free() function must not be called.      */
     numlines = (int *) mem_alloc(sizeof(ICB), "numlines");
     partition = (int *) mem_alloc(sizeof(IHBLK), "partition");
     cbval = (FLOAT *) mem_alloc(sizeof(FCB), "cbval");
     rnorm = (FLOAT *) mem_alloc(sizeof(FCB), "rnorm");
     window = (FLOAT *) mem_alloc(sizeof(FBLK), "window");
     tmn = (double *) mem_alloc(sizeof(DCB), "tmn");
     s = (FCB *) mem_alloc(sizeof(FCBCB), "s");
     lthr = (FHBLK *) mem_alloc(sizeof(F2HBLK), "lthr");
     r = (F2HBLK *) mem_alloc(sizeof(F22HBLK), "r");
     phi_sav = (F2HBLK *) mem_alloc(sizeof(F22HBLK), "phi_sav");

/*#if 0 */
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
     
     switch(lay){
	case 1: sync_flush=576; flush=384; syncsize=1024; break;
	case 2: sync_flush=480; flush=576; syncsize=1056; break;
	case 3: sync_flush=768; flush=576; syncsize=1344; break;
       default: printf("Bad lay value:(%d)",lay); exit(-1); break;
     }
/* #endif */

/* calculate HANN window coefficients */
/*   for(i=0;i<BLKSIZE;i++)  window[i]  =0.5*(1-cos(2.0*PI*i/(BLKSIZE-1.0)));*/
     for(i=0;i<BLKSIZE;i++)  window[i]  =0.5*(1-cos(2.0*PI*(i-0.5)/BLKSIZE));
     for(i=0;i<BLKSIZE_s;i++)window_s[i]=0.5*(1-cos(2.0*PI*(i-0.5)/BLKSIZE_s));
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
	for ( i=0; i<lay; i++)
  {
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
  for (j=0; j<syncsize; j++)
  {
    if (j < (sync_flush) )
      savebuf[j] = savebuf[j+flush];
    else
      savebuf[j] = *buffer++;

/**window data with HANN window***********************************************/
    if (j<BLKSIZE)
    {
      wsamp_r[j] = window[j]*((FLOAT) savebuf[j]); 
      wsamp_i[j] = 0;
    }
  }
/**Compute FFT****************************************************************/
        fft(wsamp_r,wsamp_i,energy,phi,1024, 0);
/*****************************************************************************
 * calculate the unpredictability measure, given energy[f] and phi[f]        *
 *****************************************************************************/
        for(j=0; j<HBLKSIZE; j++){
           r_prime = 2.0 * r[chn][old][j] - r[chn][oldest][j];
           phi_prime = 2.0 * phi_sav[chn][old][j] - phi_sav[chn][oldest][j];
           r[chn][new][j] = sqrt((double) energy[j]);
           phi_sav[chn][new][j] = phi[j];
	   temp1 = r[chn][new][j] * cos((double) phi[j])
		   - r_prime * cos(phi_prime);
	   temp2=r[chn][new][j] * sin((double) phi[j])
		   - r_prime * sin(phi_prime);
           temp3=r[chn][new][j] + fabs(r_prime);
           if(temp3 != 0)c[j]=sqrt(temp1*temp1+temp2*temp2)/temp3;
           else c[j] = 0;
        }
/*only update data "age" pointers after you are done with the second channel */
/*for layer 1 computations, for the layer 2 double computations, the pointers*/
/*are reset automatically on the second pass                                 */
        if(lay==2 || chn==1){
           if(new==0){new = 1; oldest = 1;}
           else {new = 0; oldest = 0;}
           if(old==0)old = 1; else old = 0;
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
/*         tb = -0.434294482*log((double) cb[j])-0.301029996; */
           tb = -0.43 *log((double) cb[j]) - 0.29 ;
           if(tb<0.0) tb=0.0; else if(tb>1.0) tb=1.0;
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
           temp1=nb[partition[j]];		 /* preliminary threshold */
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
           min_thres = fthr[j];
           sum_energy = energy[j];
           for(k=1;k<17;k++){
              if(min_thres>fthr[j+k]) min_thres = fthr[j+k];
              sum_energy += energy[j+k];
           }
           snrtmp[i][j/16] = sum_energy/(min_thres * 17.0);
           snrtmp[i][j/16] = 4.342944819 * log((double)snrtmp[i][j/16]);
        }
        for(j=208;j<(HBLKSIZE-1);j += 16){
           min_thres = 0.0;
           sum_energy = 0.0;
           for(k=0;k<17;k++){
              min_thres += fthr[j+k];
              sum_energy += energy[j+k];
           }
           snrtmp[i][j/16] = sum_energy/min_thres;
           snrtmp[i][j/16] = 4.342944819 * log((double)snrtmp[i][j/16]);
        }
/*****************************************************************************
 * End of Psychoacuostic calculation loop                                    *
 *****************************************************************************/
     }
     for(i=0; i<32; i++){
        if(lay==2)		/* if(lay==2 && chn==2) MI */
           snr32[i]=(snrtmp[0][i]>snrtmp[1][i])?snrtmp[0][i]:snrtmp[1][i];
        else snr32[i]=snrtmp[0][i];
     }
     break;

/*************************************************************************/
/**       LAYER 3                                                        */
/*************************************************************************/

   case 3:
	if ( init_L3 == 0 )
	{
	    L3para_read( sfreq,numlines,partition_l,minval,qthr_l,norm_l,s3_l,
			 partition_s,qthr_s,norm_s,SNR_s,
			 cbw_l,bu_l,bo_l,w1_l,w2_l, cbw_s,bu_s,bo_s,w1_s,w2_s );
	    init_L3 ++ ;
	}
	
	for ( j = 0; j < 21; j++ )
	    ratio_d[j] = ratio[chn][j];
	for ( j = 0; j < 12; j++ )
	    for ( i = 0; i < 3; i++ )
		ratio_ds[j][i] = ratio_s[chn][j][i];
	
	if ( chn == 0 )
	    if ( new == 0 )
	    {
		new = 1;
		old = 0;
		oldest = 1;
	    }
	    else
	    {
		new = 0;
		old = 1;
		oldest = 0;
	    }


/**********************************************************************
*  Delay signal by sync_flush=768 samples                             *
**********************************************************************/
	
	for ( j = 0; j < sync_flush; j++ ) /* for long window samples */
	    savebuf[j] = savebuf[j+flush];
	
	for ( j = sync_flush; j < syncsize; j++ )
	    savebuf[j] = *buffer++;
	
	for ( j = 0; j < BLKSIZE; j++ )
	{ /**window data with HANN window**/
	    wsamp_r[j] = window[j] * savebuf[j];  
	    wsamp_i[j] = 0.0;
	}


/**********************************************************************
*    compute unpredicatability of first six spectral lines            * 
**********************************************************************/

	fft( wsamp_r, wsamp_i, energy, phi, 1024, 1 );		/**long FFT**/
	/* energy[0--5], phi[0--5] interresting */
	for ( j = 0; j < 6; j++ )
	{	 /* calculate unpredictability measure cw */
	    r_prime = 2.0 * r[chn][old][j] - r[chn][oldest][j];
	    phi_prime = 2.0 * phi_sav[chn][old][j]-phi_sav[chn][oldest][j];
	    r[chn][new][j] = sqrt((double) energy[j]);
	    phi_sav[chn][new][j] = phi[j];
	    temp1 = r[chn][new][j] * cos((double) phi[j])
		- r_prime * cos(phi_prime);
	    temp2 = r[chn][new][j] * sin((double) phi[j])
		- r_prime * sin(phi_prime);
	    temp3=r[chn][new][j] + fabs(r_prime);
	    
	    if ( temp3 != 0.0 )
		cw[j] = sqrt( temp1*temp1+temp2*temp2 ) / temp3;
	    else
		cw[j] = 0;
	}


/**********************************************************************
*     compute unpredicatibility of next 200 spectral lines            *
**********************************************************************/ 
	for ( sblock = 0; sblock < 3; sblock++ )
	{ /**window data with HANN window**/
	    for ( j = 0, k = 128 * (2 + sblock); j < 256; j++, k++ )
	    {
		wsamp_r[j] = window_s[j]* savebuf[k]; 
		wsamp_i[j] = 0.0;
	    }							/* short FFT*/
	    
	    fft( wsamp_r, wsamp_i, &energy_s[sblock][0], &phi_s[sblock][0], 256, 2 );
        }
 
        sblock = 1;

	for ( j = 6; j < 206; j += 4 )
	{/* calculate unpredictability measure cw */
	    double r2, phi2, temp1, temp2, temp3;
	    
	    k = (j+2) >> 2;   /* k = 1 -- 51 */
	    r_prime = 2.0 * sqrt((double) energy_s[0][k])
		- sqrt((double) energy_s[2][k]);
	    phi_prime = 2.0 * phi_s[0][k] - phi_s[2][k];
	    r2 = sqrt((double) energy_s[1][k]);
	    phi2 = phi_s[1][k];
	    temp1 = r2 * cos( phi2 ) - r_prime * cos( phi_prime );
	    temp2 = r2 * sin( phi2 ) - r_prime * sin( phi_prime );
	    temp3 = r2 + fabs( r_prime );
	    if ( temp3 != 0.0 )
		cw[j] = sqrt( temp1 * temp1 + temp2 * temp2 ) / temp3;
	    else
		cw[j] = 0.0;
	    cw[j+1] = cw[j+2] = cw[j+3] = cw[j];
	}


/**********************************************************************
*    Set unpredicatiblility of remaining spectral lines to 0.4        *
**********************************************************************/
	for ( j = 206; j < HBLKSIZE; j++ )
	    cw[j] = 0.4;
	


/**********************************************************************
*    Calculate the energy and the unpredictability in the threshold   *
*    calculation partitions                                           *
**********************************************************************/

	for ( b = 0; b < CBANDS; b++ )
	{
	    eb[b] = 0.0;
	    cb[b] = 0.0;
	}
	for ( j = 0; j < HBLKSIZE; j++ )
	{
	    int tp = partition_l[j];
	    if ( tp >= 0 )
	    {
		eb[tp] += energy[j];
		cb[tp] += cw[j] * energy[j];
	    }
	}


/**********************************************************************
*      convolve the partitioned energy and unpredictability           *
*      with the spreading function, s3_l[b][k]                        *
******************************************************************** */
	
	for ( b = 0; b < CBANDS; b++ )
	{
	    ecb[b] = 0.0;
	    ctb[b] = 0.0;
	}
	if (sfreq_idx == 1) {
          sprdngf1(ecb,eb);
          sprdngf2(ctb,cb);
	} else {
	  for ( b = 0;b < CBANDS; b++ )
	  {
	    for ( k = 0; k < CBANDS; k++ )
	    {
	      if (s3_l[b][k] != 1.0) {
		ecb[b] += s3_l[b][k] * eb[k];	/* sprdngf for Layer III */
		ctb[b] += s3_l[b][k] * cb[k];
	      }
	    }
	  }
	}

	/* calculate the tonality of each threshold calculation partition */
	/* calculate the SNR in each threshhold calculation partition */

	for ( b = 0; b < CBANDS; b++ )
	{
	    double cbb,tbb;
	    if (ecb[b] != 0.0 )
                {
		cbb = ctb[b]/ecb[b];
                if (cbb <0.01) cbb = 0.01;
		cbb = log( cbb);
                }
	    else
		cbb = 0.0 ;
	    tbb = -0.299 - 0.43*cbb;  /* conv1=-0.299, conv2=-0.43 */
	    tbb = minimum( 1.0, maximum( 0.0, tbb) ) ;  /* 0<tbb<1 */
	    SNR_l[b] = maximum( minval[b], 29.0*tbb+6.0*(1.0-tbb) );
	}	/* TMN=29.0,NMT=6.0 for all calculation partitions */
	
	for ( b = 0; b < CBANDS; b++ ) /* calculate the threshold for each partition */
	    nb[b] = ecb[b] * norm_l[b] * exp( -SNR_l[b] * LN_TO_LOG10 );

	for ( b = 0; b < CBANDS; b++ )
	{ /* pre-echo control */
	    double temp_1; /* BUG of IS */
	    temp_1 = minimum( nb[b], minimum(2.0*nb_1[chn][b],16.0*nb_2[chn][b]) );
	    thr[b] = maximum( qthr_l[b], temp_1 );/* rpelev=2.0, rpelev2=16.0 */
	    nb_2[chn][b] = nb_1[chn][b];
	    nb_1[chn][b] = nb[b];
	}


	*pe = 0.0;		/*  calculate percetual entropy */
	for ( b = 0; b < CBANDS; b++ )
	{
	    double tp ;
	    tp = minimum( 0.0, log((thr[b]+1.0) / (eb[b]+1.0) ) ) ; /*not log*/
	    *pe -= numlines[b] * tp ;
	}	/* thr[b] -> thr[b]+1.0 : for non sound portition */
	
#define switch_pe  1800
        blocktype = NORM_TYPE;
	

	if ( *pe < switch_pe )
	{				/* no attack : use long blocks */
	    switch( blocktype_old[chn] ) 
	    {
	      case NORM_TYPE:
	      case STOP_TYPE:
		blocktype = NORM_TYPE;
		break;
    
	      case SHORT_TYPE:
		blocktype = STOP_TYPE;
		break;
    
	      case START_TYPE:
		fprintf( stderr, "Error in block selecting\n" );
		abort();
		break; /* problem */
	    }

	    /* threshold calculation (part 2) */
	    for ( sb = 0; sb < SBMAX_l; sb++ )
	    {
		en[sb] = w1_l[sb] * eb[bu_l[sb]] + w2_l[sb] * eb[bo_l[sb]];
		thm[sb] = w1_l[sb] *thr[bu_l[sb]] + w2_l[sb] * thr[bo_l[sb]];
		for ( b = bu_l[sb]+1; b < bo_l[sb]; b++ )
		{
		    en[sb]  += eb[b];
		    thm[sb] += thr[b];
		}
		if ( en[sb] != 0.0 )
		    ratio[chn][sb] = thm[sb]/en[sb];
		else
		    ratio[chn][sb] = 0.0;
	    }
	}
	else 
	{
	    /* attack : use short blocks */
	    blocktype = SHORT_TYPE;
	    
	    if ( blocktype_old[chn] == NORM_TYPE ) 
		blocktype_old[chn] = START_TYPE;
	    if ( blocktype_old[chn] == STOP_TYPE )
		blocktype_old[chn] = SHORT_TYPE ;
	    
	    /* threshold calculation for short blocks */
	    
	    for ( sblock = 0; sblock < 3; sblock++ )
	    {
		for ( b = 0; b < CBANDS_s; b++ )
		{
		    eb[b] = 0.0;
		    ecb[b] = 0.0;
		}
		for ( j = 0; j < HBLKSIZE_s; j++ )
		    eb[partition_s[j]] += energy_s[sblock][j];
		for ( b = 0; b < CBANDS_s; b++ )
		    for ( k = 0; k < CBANDS_s; k++ )
			ecb[b] += s3_l[b][k] * eb[k];
		for ( b = 0; b < CBANDS_s; b++ )
		{
		    nb[b] = ecb[b] * norm_l[b] * exp( (double) SNR_s[b] * LN_TO_LOG10 );
		    thr[b] = maximum (qthr_s[b],nb[b]);
		}
		for ( sb = 0; sb < SBMAX_s; sb++ )
		{
		    en[sb] = w1_s[sb] * eb[bu_s[sb]] + w2_s[sb] * eb[bo_s[sb]];
		    thm[sb] = w1_s[sb] *thr[bu_s[sb]] + w2_s[sb] * thr[bo_s[sb]];
		    for ( b = bu_s[sb]+1; b < bo_s[sb]; b++ )
		    {
			en[sb] += eb[b];
			thm[sb] += thr[b];
		    }
		    if ( en[sb] != 0.0 )
			ratio_s[chn][sb][sblock] = thm[sb]/en[sb];
		    else
			ratio_s[chn][sb][sblock] = 0.0;
		}
	    }
	} 
	
	cod_info->block_type = blocktype_old[chn];
	blocktype_old[chn] = blocktype;

	if ( cod_info->block_type == NORM_TYPE )
	    cod_info->window_switching_flag = 0;
	else
	    cod_info->window_switching_flag = 1;
	cod_info->mixed_block_flag = 0;
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
#ifdef DEBUG
#undef DEBUG
#endif


void L3para_read(double sfreq, int *numlines, int *partition_l, double *minval, double *qthr_l, double *norm_l, double (*s3_l)[63], int *partition_s, double *qthr_s, double *norm_s, double *SNR, int *cbw_l, int *bu_l, int *bo_l, double *w1_l, double *w2_l, int *cbw_s, int *bu_s, int *bo_s, double *w1_s, double *w2_s)
{
   double freq_tp;
   static double bval_l[CBANDS], bval_s[CBANDS];
   int   cbmax, cbmax_tp;
   static double s3_s[CBANDS][CBANDS];
   double *p = psy_data;
   
   char tp[256];
   int  sbmax ;
   int  i,j,k,k2,loop, part_max ;

/* Read long block data */

      for(loop=0;loop<6;loop++)
      {
	freq_tp = *p++;
	cbmax_tp = (int) *p++;
	cbmax_tp++;

	if (sfreq == freq_tp )
	  {
	     cbmax = cbmax_tp;
	     for(i=0,k2=0;i<cbmax_tp;i++)
	       {
                j = (int) *p++;
                numlines[i] = (int) *p++;
                minval[i] = *p++;
                qthr_l[i] = *p++;
                norm_l[i] = *p++;
                bval_l[i] = *p++;
	        if (j!=i)
	         { printf("please check \"psy_data\"");
		   exit(-1);
	         }
		for(k=0;k<numlines[i];k++)
		  partition_l[k2++] = i ;
		}
	   }
	   else
             p += cbmax_tp * 6;
       }

/************************************************************************
 * Now compute the spreading function, s[j][i], the value of the spread-*
 * ing function, centered at band j, for band i, store for later use    *
 ************************************************************************/
#ifdef DEBUG
	printf("freq = %f\n",sfreq);
#endif
	  part_max = cbmax ;
          for(i=0;i<part_max;i++)
	  {
	  double tempx,x,tempy,temp;
            for(j=0;j<part_max;j++)
	    {
             tempx = (bval_l[i] - bval_l[j])*1.05;
             if (j>=i) tempx = (bval_l[i] - bval_l[j])*3.0;
               else    tempx = (bval_l[i] - bval_l[j])*1.5;
/*             if (j>=i) tempx = (bval_l[j] - bval_l[i])*3.0;
               else    tempx = (bval_l[j] - bval_l[i])*1.5; */
             if(tempx>=0.5 && tempx<=2.5)
	     {
               temp = tempx - 0.5;
               x = 8.0 * (temp*temp - 2.0 * temp);
             }
             else x = 0.0;
             tempx += 0.474;
             tempy = 15.811389 + 7.5*tempx - 17.5*sqrt(1.0+tempx*tempx);
             if (tempy <= -60.0) s3_l[i][j] = 0.0;
             else                s3_l[i][j] = exp( (x + tempy)*LN_TO_LOG10 );
#ifdef DEBUG_S3
	     printf("s[%d][%d]=%f\n",i,j,s3_l[i][j]);
#endif

#ifdef DEBUGP
	     printf("j=%d i=%d tempy=%f s[i][j]=%f \n",i,j,tempy,s[i][j]);
	   minval[i] = bmax[j-1];
	   printf("minval[%d] = %f, j-1=%d %f\n",i, minval[i] , j,fthr[i]) ;
#endif
            }
          }


/* Read short block data */

      for(loop=0;loop<6;loop++)
      {
        freq_tp = *p++;
        cbmax_tp = (int) *p++;
	cbmax_tp++;

	if (sfreq == freq_tp )
	  {
	     cbmax = cbmax_tp;
	     for(i=0,k2=0;i<cbmax_tp;i++)
	       {
                j = (int) *p++;
                numlines[i] = (int) *p++;
                qthr_s[i] = *p++;
                norm_s[i] = *p++;
                SNR[i] = *p++;
                bval_s[i] = *p++;
	        if (j!=i)
	         { printf("please check \"psy_data\"");
		   exit(-1);
	         }
		for(k=0;k<numlines[i];k++)
		  partition_s[k2++] = i ;
		}
	   }
	   else
             p += cbmax_tp * 6;
       }

/************************************************************************
 * Now compute the spreading function, s[j][i], the value of the spread-*
 * ing function, centered at band j, for band i, store for later use    *
 ************************************************************************/
#ifdef DEBUG_S3
fpp=fopen("s3_s","w");
#endif
	  part_max = cbmax ;
          for(i=0;i<part_max;i++)
	  {
	  double tempx,x,tempy,temp;
            for(j=0;j<part_max;j++)
	    {
             tempx = (bval_s[i] - bval_s[j])*1.05;
             if (j>=i) tempx = (bval_s[i] - bval_s[j])*3.0;
               else    tempx = (bval_s[i] - bval_s[j])*1.5;
             if(tempx>=0.5 && tempx<=2.5)
	     {
               temp = tempx - 0.5;
               x = 8.0 * (temp*temp - 2.0 * temp);
             }
             else x = 0.0;
             tempx += 0.474;
             tempy = 15.811389 + 7.5*tempx - 17.5*sqrt(1.0+tempx*tempx);
             if (tempy <= -60.0) s3_s[i][j] = 0.0;
             else                s3_s[i][j] = exp( (x + tempy)*LN_TO_LOG10 );
#ifdef DEBUG_S3
	     fprintf(fpp,"s3_s[%d][%d]=%f\n",i,j,s3_s[i][j]);
#endif
#ifdef DEBUGP
	     printf("j=%d i=%d tempy=%f s[i][j]=%f \n",i,j,tempy,s[i][j]);
	   minval[i] = bmax[j-1];
	   printf("minval[%d] = %f, j-1=%d %f\n",i, minval[i] , j,fthr[i]) ;
#endif
            }
          }
#ifdef DEBUG_S3
	fclose(fpp);
#endif
/* Read long block data for converting threshold calculation 
   partitions to scale factor bands */

      for(loop=0;loop<6;loop++)
      {
        freq_tp = *p++;
        sbmax =  (int) *p++;
	sbmax++;

	if (sfreq == freq_tp)
	  {
	     for(i=0;i<sbmax;i++)
	      {
                j = (int) *p++;
                cbw_l[i] = (int) *p++;
                bu_l[i] = (int) *p++;
                bo_l[i] = (int) *p++;
                w1_l[i] = (double) *p++;
                w2_l[i] = (double) *p++;
	        if (j!=i)
	         { printf("30:please check \"psy_data\"\n");
		   exit(-1);
	         }
	        if (i!=0)
		 if ( (bo_l[i] != (bu_l[i]+cbw_l[i])) ||
				 (fabs(1.0-w1_l[i]-w2_l[i-1]) > 0.01 ) )
	         { printf("31:please check \"psy_data.\"\n");
		   exit(-1);
	         }
	      }
	   }
	   else
             p += sbmax * 6;
       }

/* Read short block data for converting threshold calculation 
   partitions to scale factor bands */

      for(loop=0;loop<6;loop++)
      {
        freq_tp = *p++;
        sbmax = (int) *p++;
	sbmax++;

	if (sfreq == freq_tp)
	  {
	     for(i=0;i<sbmax;i++)
	      {
                j = (int) *p++;
                cbw_s[i] = (int) *p++;
                bu_s[i] = (int) *p++;
                bo_s[i] = (int) *p++;
                w1_s[i] = *p++;
                w2_s[i] = *p++;
	        if (j!=i)
	         { printf("30:please check \"psy_data\"\n");
		   exit(-1);
	         }
	        if (i!=0)
		 if ( (bo_s[i] != (bu_s[i]+cbw_s[i])) ||
				 (fabs(1.0-w1_s[i]-w2_s[i-1]) > 0.01 ) )
	         { printf("31:please check \"psy_data.\"\n");
		   exit(-1);
	         }
	      }
	   }
	   else
	     p += sbmax * 6;
       }

}

 static int s3ind[CBANDS][2] = {
   {0,2},
   {0,3},
   {0,4},
   {0,5},
   {0,6},
   {0,7},
   {0,8},
   {0,9},
   {0,10},
   {0,11},
   {0,12},
   {1,14},
   {1,14},
   {2,15},
   {3,15},
   {5,16},
   {6,17},
   {7,19},
   {9,20},
   {10,21},
   {11,22},
   {12,23},
   {14,24},
   {15,25},
   {15,27},
   {16,28},
   {16,28},
   {17,29},
   {18,30},
   {19,31},
   {19,32},
   {20,34},
   {21,35},
   {22,36},
   {22,36},
   {23,37},
   {24,38},
   {25,39},
   {26,41},
   {27,42},
   {28,43},
   {29,44},
   {30,45},
   {31,46},
   {32,47},
   {33,48},
   {34,49},
   {35,50},
   {36,51},
   {37,52},
   {37,53},
   {38,54},
   {39,55},
   {40,56},
   {41,57},
   {42,58},
   {43,59},
   {44,60},
   {45,61},
   {46,62},
   {47,62},
   {48,62},
   {48,62},
 };

void sprdngf1(dest,source)
FLOAT *dest;
double *source;
{
static int mfcinit = 0;
int b,k;
 
        for ( b = 0;b < CBANDS; b++ )
            for ( k = s3ind[b][0]; k <= s3ind[b][1]; k++ )
                dest[b] += s3_l[b][k] * source[k];      
}

void sprdngf2(dest,source)
double *dest;
float *source;
{
static int mfcinit = 0;
int b,k;

        for ( b = 0;b < CBANDS; b++ )
            for ( k = s3ind[b][0]; k <= s3ind[b][1]; k++ )
                dest[b] += s3_l[b][k] * source[k];      
}
