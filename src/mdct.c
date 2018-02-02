/**********************************************************************
 * ISO MPEG Audio Subgroup Software Simulation Group (1996)
 * ISO 13818-3 MPEG-2 Audio Encoder - Lower Sampling Frequency Extension
 *
 * $Id: mdct.c,v 1.1 1996/02/14 04:04:23 rowlands Exp $
 *
 * $Log: mdct.c,v $
 * Revision 1.1  1996/02/14 04:04:23  rowlands
 * Initial revision
 *
 * Received from Mike Coleman
 **********************************************************************/

#include "common.h"
#include "l3side.h"
#include "mdct.h"

double ca[8], cs[8];

/*
  This is table B.9: coefficients for aliasing reduction
*/
static double c[8] = { -0.6,-0.535,-0.33,-0.185,-0.095,-0.041,-0.0142, -0.0037 };

void mdct_sub( L3SBS (*sb_sample), double (*mdct_freq)[2][576], int stereo, III_side_info_t *l3_side, int mode_gr )
{
    gr_info *cod_info;
    double mdct_in[36];
    int ch,gr,band,k,j;
    double bu,bd;
    static int init = 0;
    int	block_type;
    double (*mdct_enc)[2][32][18] = (double (*)[2][32][18]) mdct_freq;
    
    if ( init == 0 )
    {
	/* prepare the aliasing reduction butterflies */
	for ( k = 0; k < 8; k++ )
	{
	    double sq;
	    sq = sqrt( 1.0 + c[k] * c[k] );
	    ca[k] = c[k] / sq;
	    cs[k] = 1.0 / sq;
	}
	init++;
    }
    
    for ( gr = 0; gr < mode_gr; gr++ )
	for ( ch = 0; ch < stereo; ch++ )
	{
	    cod_info = (gr_info *) &(l3_side->gr[gr].ch[ch]) ;
	    block_type = cod_info->block_type;
	    
	    /*
	      Compensate for inversion in the analysis filter
	    */
	    for ( band = 0; band < 32; band++ )
		for ( k = 0; k < 18; k++ )
		    if ( (band & 1) && (k & 1) )
			(*sb_sample)[ch][gr+1][k][band] *= -1.0;
	    
	    /*
	      Perform imdct of 18 previous subband samples
	      + 18 current subband samples
	    */
	    for ( band = 0; band < 32; band++ )
	    {
		for ( k = 0; k < 18; k++ )
		{
		    mdct_in[k]    = (*sb_sample)[ch][ gr ][k][band];
		    mdct_in[k+18] = (*sb_sample)[ch][gr+1][k][band];
		}
		if ( cod_info->mixed_block_flag && (band < 2) )
		    block_type = 0;
		
		mdct( mdct_in, &mdct_enc[gr][ch][band][0], block_type );
	    }
	    
	    /*
	      Perform aliasing reduction butterfly
	      on long blocks
	    */
	    if ( block_type != 2 )
		for ( band = 0; band < 31; band++ )
		    for ( k = 0; k < 8; k++ )
		    {
			bu = mdct_enc[gr][ch][band][17-k] * cs[k] + mdct_enc[gr][ch][band+1][k] * ca[k];
			bd = mdct_enc[gr][ch][band+1][k] * cs[k] - mdct_enc[gr][ch][band][17-k] * ca[k];
			mdct_enc[gr][ch][band][17-k] = bu;
			mdct_enc[gr][ch][band+1][k]  = bd;
		    }
	    
	}
    
    /*
      Save latest granule's subband samples to be used in
      the next mdct call
    */
    for ( ch = 0; ch < stereo; ch++ )
	for ( j = 0; j < 18; j++ )
	    for ( band = 0; band < 32; band++ )
		(*sb_sample)[ch][0][j][band] = (*sb_sample)[ch][mode_gr][j][band];
}

void mdct( double *in, double *out, int block_type )
{
/*-------------------------------------------------------------------*/
/*                                                                   */
/*   Function: Calculation of the MDCT                               */
/*   In the case of long blocks ( block_type 0,1,3 ) there are       */
/*   36 coefficents in the time domain and 18 in the frequency       */
/*   domain.                                                         */
/*   In the case of short blocks (block_type 2 ) there are 3         */
/*   transformations with short length. This leads to 12 coefficents */
/*   in the time and 6 in the frequency domain. In this case the     */
/*   results are stored side by side in the vector out[].            */
/*                                                                   */
/*   New layer3                                                      */
/*                                                                   */
/*-------------------------------------------------------------------*/

  int l,k,i,m,N;
  double sum;
  static double win[4][36];
  static int init = 0;
  static double cos_s[6][12], cos_l[18][36];
  static double fin[36];

  if ( init == 0 )
  {
    /* type 0 */
    for ( i = 0; i < 36; i++ )
      win[0][i] = sin( PI/36 * (i + 0.5) );
    /* type 1*/
    for ( i = 0; i < 18; i++ ) 
      win[1][i] = sin( PI/36 * (i + 0.5) );
    for ( i = 18; i < 24; i++ )
      win[1][i] = 1.0;
    for ( i = 24; i < 30; i++ )
      win[1][i] = sin( PI/12 * ( i + 0.5 - 18) );
    for ( i = 30; i < 36; i++ )
      win[1][i] = 0.0;
    /* type 3*/
    for ( i = 0; i < 6; i++ )
      win[3][i] = 0.0;
    for ( i = 6; i < 12; i++ ) 
      win[3][i] = sin( PI/12 * (i + 0.5 - 6) );
    for ( i = 12; i < 18; i++ )
      win[3][i] = 1.0;
    for ( i = 18; i < 36; i++ )
      win[3][i] = sin( PI/36 * (i + 0.5) );
    /* type 2*/
    for ( i = 0; i < 12; i++ )
    win[2][i] = sin( PI/12 * (i + 0.5) );
    for ( i = 12; i < 36; i++ )
      win[2][i] = 0.0;

    N = 12;
    for ( m = 0; m < N / 2; m++ )
      for ( k = 0; k < N; k++ )
        cos_s[m][k] = cos( (PI /(2 * N)) * (2 * k + 1 + N / 2) *
                     (2 * m + 1) ) / (N / 4);

    N = 36;
    for ( m = 0; m < N / 2; m++ )
      for ( k = 0; k < N; k++ )
        cos_l[m][k] = cos( (PI / (2 * N)) * (2 * k + 1 + N / 2) *
                     (2 * m + 1) ) / (N / 4);

    init++;
  }

  if ( block_type == 2 )
  {
    N = 12;
    for ( l = 0; l < 3; l++ )
    {
      for ( m = 0; m < N / 2; m++ )
      {
        for ( sum = 0.0, k = 0; k < N; k++ )
          sum += win[block_type][k] * in[k + 6 * l + 6] * cos_s[m][k];
        out[ 3 * m + l] = sum;
      }
    }
  }
  else
  {
  if (block_type!=0) { /* then it's a short->long or long->short block */
    N = 36;
    for (k=0;k<N;k++)
      fin[k]=win[block_type][k]*in[k];
    for ( m = 0; m < N / 2; m++ )
    {
      for ( sum = 0.0, k = 0; k < N; k++ )
        sum += fin[k] * cos_l[m][k];
      out[m] = sum;
    }
  }
	else /* block_type is 0.  This is a long->long transform */
    {
      N=36;
      for (k=0;k<N;k++)
	fin[k]=win[0][k]*in[k];   /* do the f(k)*in(k) first. and save the results */

      /* 0 */  
      sum = ( fin[0]-fin[17] ) * cos_l[0][0]; /* 17 */
      sum += ( fin[1]-fin[16] ) * cos_l[0][1]; /* 15 */
      sum += ( fin[2]-fin[15] ) * cos_l[0][2]; /* 13 */
      sum += ( fin[3]-fin[14] ) * cos_l[0][3]; /* 11 */
      sum += ( fin[4]-fin[13] ) * cos_l[0][4]; /* 9 */
      sum += ( fin[5]-fin[12] ) * cos_l[0][5]; /* 7  */
      sum += ( fin[6]-fin[11] ) * cos_l[0][6]; /* 5  */
      sum += ( fin[7]-fin[10] ) * cos_l[0][7]; /* 3 */
      sum += ( fin[8]-fin[9] ) * cos_l[0][8]; /* 1  */
      sum += ( -fin[18]-fin[35] ) * -cos_l[0][18]; /* 19*/
      sum += ( -fin[19]-fin[34] ) * -cos_l[0][19]; /* 21 */
      sum += ( -fin[20]-fin[33] ) * -cos_l[0][20]; /* 23 */
      sum += ( -fin[21]-fin[32] ) * -cos_l[0][21]; /*25 */
      sum += ( -fin[22]-fin[31] ) * -cos_l[0][22]; /* 27 */
      sum += ( -fin[23]-fin[30] ) * -cos_l[0][23]; /* 29 */
      sum += ( -fin[24]-fin[29] ) * -cos_l[0][24]; /* 31*/
      sum += ( -fin[25]-fin[28] ) * -cos_l[0][25]; /* 33 */
      sum += ( -fin[26]-fin[27] ) * -cos_l[0][26]; /* 35 */    
      out[0]=sum; 


      /* 1 */
      sum = ( -fin[2]-fin[3]+fin[14]+fin[15]-fin[26]-fin[27] ) * -cos_l[1][2];  /* mfc=3 2*/
      sum += ( -fin[1]-fin[4]+fin[13]+fin[16]-fin[25]-fin[28] ) * -cos_l[1][1];  /* mfc=9 1*/
      sum += ( -fin[0]-fin[5]+fin[12]+fin[17]-fin[24]-fin[29] ) * -cos_l[1][0];  /* mfc=15 0*/
      sum += ( -fin[6]+fin[11]+fin[18]-fin[23]-fin[30]+fin[35] ) * -cos_l[1][6]; /* mfc=21 6*/
      sum += ( -fin[7]+fin[10]+fin[19]-fin[22]-fin[31]+fin[34] ) * -cos_l[1][7]; /* mfc=27 7*/
      sum += ( -fin[8]+fin[9]+fin[20]-fin[21]-fin[32]+fin[33] ) * -cos_l[1][8]; /* mfc = 28 8*/
      out[1]=sum;

      /* 2 */
      sum = ( -fin[0]+fin[17] ) * -cos_l[2][0];  /* mfc 23 */
      sum += ( -fin[1]+fin[16] ) * -cos_l[2][1];  /* mfc 33 */
      sum += ( fin[2]-fin[15] ) * cos_l[2][2];   /* mfc 29 */
      sum += ( fin[3]-fin[14] ) * cos_l[2][3];   /* mfc 19 */
      sum += ( fin[4]-fin[13] ) * cos_l[2][4];   /* mfc 9  */
      sum += ( fin[5]-fin[12] ) * cos_l[2][5];   /* mfc 1  */
      sum += ( fin[6]-fin[11] ) * cos_l[2][6];   /* mfc 11 */
      sum += ( fin[7]-fin[10] ) * cos_l[2][7];   /* mfc 21 */
      sum += ( fin[8]-fin[9] ) * cos_l[2][8];    /* mfc 31 */
      sum += ( fin[18]+fin[35] ) * cos_l[2][18]; /* mfc 13 */
      sum += ( fin[19]+fin[34] ) * cos_l[2][19]; /* mfc 3  */
      sum += ( fin[20]+fin[33] ) * cos_l[2][20]; /* mfc 7  */
      sum += ( fin[21]+fin[32] ) * cos_l[2][21]; /* mfc 17 */
      sum += ( fin[22]+fin[31] ) * cos_l[2][22]; /* mfc 27 */
      sum += ( -fin[23]-fin[30] ) * -cos_l[2][23]; /* mfc 35 */
      sum += ( -fin[24]-fin[29] ) * -cos_l[2][24]; /* mfc 25 */
      sum += ( -fin[25]-fin[28] ) * -cos_l[2][25]; /* mfc 15 */
      sum += ( -fin[26]-fin[27] ) * -cos_l[2][26]; /* mfc 5 */
      out[2]=sum;

      /* 3 */
      sum = ( fin[0]-fin[17] ) * cos_l[3][0]; /* mfc 11 */
      sum += ( fin[1]-fin[16] ) * cos_l[3][1]; /* mfc 3 */
      sum += ( fin[2]-fin[15] ) * cos_l[3][2]; /* 17 */
      sum += ( fin[3]-fin[14] ) * cos_l[3][3]; /* 31 */
      sum += ( -fin[4]+fin[13] ) * -cos_l[3][4]; /* 27 */
      sum += ( -fin[5]+fin[12] ) * -cos_l[3][5]; /* 13 */
      sum += ( -fin[6]+fin[11] ) * -cos_l[3][6]; /* 1 */
      sum += ( -fin[7]+fin[10] ) * -cos_l[3][7]; /* 15 */
      sum += ( -fin[8]+fin[9] ) * -cos_l[3][8]; /* 29 */
      sum += ( -fin[18]-fin[35] ) * -cos_l[3][18]; /* 25 */
      sum += ( fin[19]+fin[34] ) * cos_l[3][19]; /* 33 */
      sum += ( fin[20]+fin[33] ) * cos_l[3][20]; /* 19 */
      sum += ( fin[21]+fin[32] ) * cos_l[3][21]; /* 5 */
      sum += ( fin[22]+fin[31] ) * cos_l[3][22]; /* 9 */
      sum += ( fin[23]+fin[30] ) * cos_l[3][23]; /* 23 */
      sum += ( -fin[24]-fin[29] ) * -cos_l[3][24]; /* 35 */
      sum += ( -fin[25]-fin[28] ) * -cos_l[3][25]; /* 21 */
      sum += ( -fin[26]-fin[27] ) * -cos_l[3][26]; /* 7 */
      out[3]=sum;

      /* 4 */
      /* mfc 27 */
      sum = ( fin[0]-fin[1]-fin[4]+fin[5]+fin[8]-fin[9]-fin[12]+fin[13]+fin[16]-fin[17] \
	       -fin[20]+fin[21]+fin[24]-fin[25]-fin[28]+fin[29]+fin[32]-fin[33] ) * cos_l[4][0];
      /* mfc 9 */
      sum += ( -fin[2]-fin[3]+fin[6]+fin[7]-fin[10]-fin[11]+fin[14]+fin[15]-fin[18]-fin[19] \
	       +fin[22]+fin[23]-fin[26]-fin[27]+fin[30]+fin[31]-fin[34]-fin[35] ) * cos_l[4][6];
      out[4]=sum;

      /* 5 */
      sum = ( -fin[0]+fin[17] ) * -cos_l[5][0]; /* 7 */
      sum += ( -fin[1]+fin[16] ) * -cos_l[5][1]; /* 15 */
      sum += ( fin[2]-fin[15] ) * cos_l[5][2]; /* 35 */
      sum += ( fin[3]-fin[14] ) * cos_l[5][3]; /* 13 */
      sum += ( fin[4]-fin[13] ) * cos_l[5][4]; /* 9 */
      sum += ( fin[5]-fin[12] ) * cos_l[5][5]; /* 31 */
      sum += ( -fin[6]+fin[11] ) * -cos_l[5][6]; /* 19 */
      sum += ( -fin[7]+fin[10] ) * -cos_l[5][7]; /* 3 */
      sum += ( -fin[8]+fin[9] ) * -cos_l[5][8]; /* 8 */
      sum += ( fin[18]+fin[35] ) * cos_l[5][18]; /* 29 */
      sum += ( -fin[19]-fin[34] ) * -cos_l[5][19]; /* 21 */
      sum += ( -fin[20]-fin[33] ) * -cos_l[5][20]; /* 1 */
      sum += ( -fin[21]-fin[32] ) * -cos_l[5][21]; /* 23 */
      sum += ( fin[22]+fin[31] ) * cos_l[5][22]; /* 27 */
      sum += ( fin[23]+fin[30] ) * cos_l[5][23]; /* 5 */
      sum += ( fin[24]+fin[29] ) * cos_l[5][24]; /* 17 */
      sum += ( -fin[25]-fin[28] ) * -cos_l[5][25]; /* 33 */
      sum += ( -fin[26]-fin[27] ) * -cos_l[5][26]; /* 11 */
      out[5]=sum;

      /* 6 */
      sum = ( -fin[0]+fin[17] ) * -cos_l[6][0]; /* 31 */
      sum += ( fin[1]-fin[16] ) * cos_l[6][1]; /* 15 */
      sum += ( fin[2]-fin[15] ) * cos_l[6][2]; /* 11 */
      sum += ( -fin[3]+fin[14] ) * -cos_l[6][3]; /* 35 */
      sum += ( -fin[4]+fin[13] ) * -cos_l[6][4]; /* 9 */
      sum += ( -fin[5]+fin[12] ) * -cos_l[6][5]; /* 17 */
      sum += ( fin[6]-fin[11] ) * cos_l[6][6]; /* 29 */
      sum += ( fin[7]-fin[10] ) * cos_l[6][7]; /* 3 */
      sum += ( fin[8]-fin[9] ) * cos_l[6][8]; /* 23 */
      sum += ( fin[18]+fin[35] ) * cos_l[6][18]; /* 5 */
      sum += ( fin[19]+fin[34] ) * cos_l[6][19]; /* 21 */
      sum += ( -fin[20]-fin[33] ) * -cos_l[6][20]; /* 25 */
      sum += ( -fin[21]-fin[32] ) * -cos_l[6][21]; /* 1 */
      sum += ( -fin[22]-fin[31] ) * -cos_l[6][22]; /* 27 */
      sum += ( fin[23]+fin[30] ) * cos_l[6][23]; /* 19 */
      sum += ( fin[24]+fin[29] ) * cos_l[6][24]; /* 7 */
      sum += ( fin[25]+fin[28] ) * cos_l[6][25]; /* 33 */
      sum += ( -fin[26]-fin[27] ) * -cos_l[6][26]; /* 13 */    
      out[6]=sum;

      /* 7 */
      sum = ( fin[0]+fin[5]-fin[12]-fin[17]+fin[24]+fin[29] ) * cos_l[7][0]; /* 3 */
      sum += ( fin[1]+fin[4]-fin[13]-fin[16]+fin[25]+fin[28] ) * cos_l[7][1]; /* 27*/
      sum += ( -fin[2]-fin[3]+fin[14]+fin[15]-fin[26]-fin[27] ) * -cos_l[7][2]; /* 15 */
      sum += ( fin[6]-fin[11]-fin[18]+fin[23]+fin[30]-fin[35] ) * cos_l[7][6]; /* 33 */
      sum += ( -fin[7]+fin[10]+fin[19]-fin[22]-fin[31]+fin[34] ) * -cos_l[7][7]; /* 9 */
      sum += ( -fin[8]+fin[9]+fin[20]-fin[21]-fin[32]+fin[33] ) * -cos_l[7][8]; /* 21 */
      out[7]=sum;

      /* 8 */
      sum = (  fin[0]-fin[17] ) * cos_l[8][0]; /* 35 */
      sum += ( -fin[1]+fin[16] ) * -cos_l[8][1]; /* 3  */
      sum += ( -fin[2]+fin[15] ) * -cos_l[8][2]; /* 31 */
      sum += ( fin[3]-fin[14] ) * cos_l[8][3]; /*  7 */
      sum += ( fin[4]-fin[13] ) * cos_l[8][4]; /* 27*/
      sum += ( -fin[5]+fin[12] ) * -cos_l[8][5]; /* 11 */
      sum += ( -fin[6]+fin[11] ) * -cos_l[8][6]; /* 23 */
      sum += ( fin[7]-fin[10] ) * cos_l[8][7]; /* 15*/
      sum += ( fin[8]-fin[9] ) * cos_l[8][8]; /* 19 */
      sum += ( -fin[18]-fin[35] ) * -cos_l[8][18]; /* 1 */
      sum += ( -fin[19]-fin[34] ) * -cos_l[8][19]; /* 33 */
      sum += ( fin[20]+fin[33] ) * cos_l[8][20]; /*  5 */
      sum += ( fin[21]+fin[32] ) * cos_l[8][21]; /* 29*/
      sum += ( -fin[22]-fin[31] ) * -cos_l[8][22]; /*  9 */
      sum += ( -fin[23]-fin[30] ) * -cos_l[8][23]; /* 25 */
      sum += ( fin[24]+fin[29] ) * cos_l[8][24]; /* 13*/
      sum += ( fin[25]+fin[28] ) * cos_l[8][25]; /* 21 */
      sum += ( -fin[26]-fin[27] ) * -cos_l[8][26]; /* 17 */    
      out[8]=sum;

      /* 9 */
      sum = ( -fin[0]+fin[17] ) * -cos_l[9][0]; /* 1  */
      sum += ( fin[1]-fin[16] ) * cos_l[9][1]; /* 33 */
      sum += ( fin[2]-fin[15] ) * cos_l[9][2]; /*  5 */
      sum += ( -fin[3]+fin[14] ) * -cos_l[9][3]; /*  29*/
      sum += ( -fin[4]+fin[13] ) * -cos_l[9][4]; /* 9 */
      sum += ( fin[5]-fin[12] ) * cos_l[9][5]; /* 25 */
      sum += ( fin[6]-fin[11] ) * cos_l[9][6]; /* 13 */
      sum += ( -fin[7]+fin[10] ) * -cos_l[9][7]; /* 21*/
      sum += ( -fin[8]+fin[9] ) * -cos_l[9][8]; /* 17 */
      sum += ( -fin[18]-fin[35] ) * -cos_l[9][18]; /* 35*/
      sum += ( -fin[19]-fin[34] ) * -cos_l[9][19]; /* 3  */
      sum += ( fin[20]+fin[33] ) * cos_l[9][20]; /* 31 */
      sum += ( fin[21]+fin[32] ) * cos_l[9][21]; /* 7 */
      sum += ( -fin[22]-fin[31] ) * -cos_l[9][22]; /* 27 */
      sum += ( -fin[23]-fin[30] ) * -cos_l[9][23]; /* 11 */
      sum += ( fin[24]+fin[29] ) * cos_l[9][24]; /* 23*/
      sum += ( fin[25]+fin[28] ) * cos_l[9][25]; /* 15 */
      sum += ( -fin[26]-fin[27] ) * -cos_l[9][26]; /* 19 */    
      out[9]=sum; 

      /* 10 */
      sum = ( fin[0]+fin[5]-fin[12]-fin[17]+fin[24]+fin[29] ) * cos_l[10][0]; /* 33 */
      sum += ( fin[1]+fin[4]-fin[13]-fin[16]+fin[25]+fin[28] ) * cos_l[10][1]; /* 9 */
      sum += ( -fin[2]-fin[3]+fin[14]+fin[15]-fin[26]-fin[27] ) * -cos_l[10][2]; /* 21 */
      sum += ( -fin[6]+fin[11]+fin[18]-fin[23]-fin[30]+fin[35] ) * -cos_l[10][6]; /* 3 */
      sum += ( fin[7]-fin[10]-fin[19]+fin[22]+fin[31]-fin[34] ) * cos_l[10][7]; /* 27 */
      sum += ( fin[8]-fin[9]-fin[20]+fin[21]+fin[32]-fin[33] ) * cos_l[10][8]; /* 15 */
      out[10]=sum;

      /* 11 */
      sum = ( fin[0]-fin[17] ) * cos_l[11][0]; /* 5  */
      sum += ( -fin[1]+fin[16] ) * -cos_l[11][1]; /* 21 */
      sum += ( -fin[2]+fin[15] ) * -cos_l[11][2]; /* 25 */
      sum += ( fin[3]-fin[14] ) * cos_l[11][3]; /*  1 */
      sum += ( -fin[4]+fin[13] ) * -cos_l[11][4]; /*27 */
      sum += ( -fin[5]+fin[12] ) * -cos_l[11][5]; /* 19 */
      sum += ( fin[6]-fin[11] ) * cos_l[11][6]; /* 7  */
      sum += ( fin[7]-fin[10] ) * cos_l[11][7]; /* 33*/
      sum += ( -fin[8]+fin[9] ) * -cos_l[11][8]; /* 13 */
      sum += ( fin[18]+fin[35] ) * cos_l[11][18]; /* 31*/
      sum += ( fin[19]+fin[34] ) * cos_l[11][19]; /* 15 */
      sum += ( -fin[20]-fin[33] ) * -cos_l[11][20]; /* 11 */
      sum += ( -fin[21]-fin[32] ) * -cos_l[11][21]; /*35 */
      sum += ( fin[22]+fin[31] ) * cos_l[11][22]; /*  9 */
      sum += ( -fin[23]-fin[30] ) * -cos_l[11][23]; /* 17 */
      sum += ( -fin[24]-fin[29] ) * -cos_l[11][24]; /* 29*/
      sum += ( fin[25]+fin[28] ) * cos_l[11][25]; /* 3  */
      sum += ( -fin[26]-fin[27] ) * -cos_l[11][26]; /* 23 */    
      out[11]=sum; 

      /* 12 */
      sum = ( -fin[0]+fin[17] ) * -cos_l[12][0]; /* 29 */
      sum += ( -fin[1]+fin[16] ) * -cos_l[12][1]; /* 21 */
      sum += ( fin[2]-fin[15] ) * cos_l[12][2]; /* 1  */
      sum += ( -fin[3]+fin[14] ) * -cos_l[12][3]; /* 23 */
      sum += ( -fin[4]+fin[13] ) * -cos_l[12][4]; /*27 */
      sum += ( fin[5]-fin[12] ) * cos_l[12][5]; /* 5  */
      sum += ( -fin[6]+fin[11] ) * -cos_l[12][6]; /* 17 */
      sum += ( -fin[7]+fin[10] ) * -cos_l[12][7]; /* 33*/
      sum += ( fin[8]-fin[9] ) * cos_l[12][8]; /* 11 */
      sum += ( -fin[18]-fin[35] ) * -cos_l[12][18]; /* 7 */
      sum += ( fin[19]+fin[34] ) * cos_l[12][19]; /* 15 */
      sum += ( fin[20]+fin[33] ) * cos_l[12][20]; /* 35 */
      sum += ( -fin[21]-fin[32] ) * -cos_l[12][21]; /*13 */
      sum += ( fin[22]+fin[31] ) * cos_l[12][22]; /*  9 */
      sum += ( -fin[23]-fin[30] ) * -cos_l[12][23]; /* 31 */
      sum += ( -fin[24]-fin[29] ) * -cos_l[12][24]; /* 19*/
      sum += ( fin[25]+fin[28] ) * cos_l[12][25]; /* 3  */
      sum += ( -fin[26]-fin[27] ) * -cos_l[12][26]; /* 25 */    
      out[12]=sum; 

      /* 13 */
      sum = ( -fin[0]+fin[1]+fin[4]-fin[5]-fin[8]+fin[9]+fin[12]-fin[13]-fin[16]+fin[17]+fin[20]-fin[21] \
	      -fin[24]+fin[25]+fin[28]-fin[29]-fin[32]+fin[33] ) * -cos_l[13][0]; /* 9 */
      sum += ( -fin[2]-fin[3]+fin[6]+fin[7]-fin[10]-fin[11]+fin[14]+fin[15]-fin[18]-fin[19]+fin[22] \
	       +fin[23]-fin[26]-fin[27]+fin[30]+fin[31]-fin[34]-fin[35] ) * -cos_l[13][2]; /* 27 */
      out[13]=sum;

      /* 14 */
      sum = ( fin[0]-fin[17] ) * cos_l[14][0]; /* 25 */
      sum += ( fin[1]-fin[16] ) * cos_l[14][1]; /* 33 */
      sum += ( -fin[2]+fin[15] ) * -cos_l[14][2]; /* 19 */
      sum += ( fin[3]-fin[14] ) * cos_l[14][3]; /* 5  */
      sum += ( -fin[4]+fin[13] ) * -cos_l[14][4]; /* 9 */
      sum += ( fin[5]-fin[12] ) * cos_l[14][5]; /* 23 */
      sum += ( fin[6]-fin[11] ) * cos_l[14][6]; /* 35 */
      sum += ( -fin[7]+fin[10] ) * -cos_l[14][7]; /* 21*/
      sum += ( fin[8]-fin[9] ) * cos_l[14][8]; /* 7  */
      sum += ( fin[18]+fin[35] ) * cos_l[14][18]; /* 11*/
      sum += ( -fin[19]-fin[34] ) * -cos_l[14][19]; /* 3  */
      sum += ( fin[20]+fin[33] ) * cos_l[14][20]; /* 17 */
      sum += ( -fin[21]-fin[32] ) * -cos_l[14][21]; /*31 */
      sum += ( -fin[22]-fin[31] ) * -cos_l[14][22]; /* 27 */
      sum += ( fin[23]+fin[30] ) * cos_l[14][23]; /* 13 */
      sum += ( -fin[24]-fin[29] ) * -cos_l[14][24]; /* 1 */
      sum += ( fin[25]+fin[28] ) * cos_l[14][25]; /* 15 */
      sum += ( -fin[26]-fin[27] ) * -cos_l[14][26]; /* 29 */    
      out[14]=sum; 

      /* 15 */
      sum = ( fin[0]-fin[17] ) * cos_l[15][0]; /* 13 */
      sum += ( -fin[1]+fin[16] ) * -cos_l[15][1]; /* 3  */
      sum += ( fin[2]-fin[15] ) * cos_l[15][2]; /* 7  */
      sum += ( -fin[3]+fin[14] ) * -cos_l[15][3]; /* 17 */
      sum += ( fin[4]-fin[13] ) * cos_l[15][4]; /* 27*/
      sum += ( fin[5]-fin[12] ) * cos_l[15][5]; /* 35 */
      sum += ( -fin[6]+fin[11] ) * -cos_l[15][6]; /* 25 */
      sum += ( fin[7]-fin[10] ) * cos_l[15][7]; /* 15*/
      sum += ( -fin[8]+fin[9] ) * -cos_l[15][8]; /* 5  */
      sum += ( fin[18]+fin[35] ) * cos_l[15][18]; /* 23*/
      sum += ( -fin[19]-fin[34] ) * -cos_l[15][19]; /* 33 */
      sum += ( -fin[20]-fin[33] ) * -cos_l[15][20]; /* 29 */
      sum += ( fin[21]+fin[32] ) * cos_l[15][21]; /*19 */
      sum += ( -fin[22]-fin[31] ) * -cos_l[15][22]; /* 9  */
      sum += ( fin[23]+fin[30] ) * cos_l[15][23]; /* 1  */
      sum += ( -fin[24]-fin[29] ) * -cos_l[15][24]; /* 11*/
      sum += ( fin[25]+fin[28] ) * cos_l[15][25]; /* 21 */
      sum += ( -fin[26]-fin[27] ) * -cos_l[15][26]; /* 31 */    
      out[15]=sum; 

      /* 16 */
      sum = ( -fin[0]-fin[5]+fin[12]+fin[17]-fin[24]-fin[29] ) * -cos_l[16][0]; /* 21 */
      sum += ( fin[1]+fin[4]-fin[13]-fin[16]+fin[25]+fin[28] ) * cos_l[16][1];  /* 27 */
      sum += ( -fin[2]-fin[3]+fin[14]+fin[15]-fin[26]-fin[27] ) * -cos_l[16][2]; /* 33 */
      sum += ( fin[6]-fin[11]-fin[18]+fin[23]+fin[30]-fin[35] ) * cos_l[16][6]; /* 15 */
      sum += ( -fin[7]+fin[10]+fin[19]-fin[22]-fin[31]+fin[34] ) * -cos_l[16][7]; /* 9 */
      sum += ( fin[8]-fin[9]-fin[20]+fin[21]+fin[32]-fin[33] ) * cos_l[16][8]; /* 3 */
      out[16]=sum;

      /* 17 */
      sum = ( -fin[0]+fin[17] ) * -cos_l[17][0]; /* 17 */
      sum += ( fin[1]-fin[16] ) * cos_l[17][1]; /* 15 */
      sum += ( -fin[2]+fin[15] ) * -cos_l[17][2]; /* 13 */
      sum += ( fin[3]-fin[14] ) * cos_l[17][3]; /* 11 */
      sum += ( -fin[4]+fin[13] ) * -cos_l[17][4]; /* 9 */
      sum += ( fin[5]-fin[12] ) * cos_l[17][5]; /* 7  */
      sum += ( -fin[6]+fin[11] ) * -cos_l[17][6]; /* 5  */
      sum += ( fin[7]-fin[10] ) * cos_l[17][7]; /* 3 */
      sum += ( -fin[8]+fin[9] ) * -cos_l[17][8]; /* 1  */
      sum += ( -fin[18]-fin[35] ) * -cos_l[17][18]; /* 19*/
      sum += ( fin[19]+fin[34] ) * cos_l[17][19]; /* 21 */
      sum += ( -fin[20]-fin[33] ) * -cos_l[17][20]; /* 23 */
      sum += ( fin[21]+fin[32] ) * cos_l[17][21]; /*25 */
      sum += ( -fin[22]-fin[31] ) * -cos_l[17][22]; /* 27 */
      sum += ( fin[23]+fin[30] ) * cos_l[17][23]; /* 29 */
      sum += ( -fin[24]-fin[29] ) * -cos_l[17][24]; /* 31*/
      sum += ( fin[25]+fin[28] ) * cos_l[17][25]; /* 33 */
      sum += ( -fin[26]-fin[27] ) * -cos_l[17][26]; /* 35 */    
      out[17]=sum; 
    }
  }
}

void
delay( double (*xr)[2][576], int stereo )
{
    static double xr_buff[2][576];
    double xr_buff2[2][576];
    unsigned int i,j;
    
    for (i=0;i<stereo;i++)
    {
	for (j=0;j<576;j++) xr_buff2[i][j] = xr_buff[i][j];
	for (j=0;j<576;j++) xr_buff[i][j]  = xr[1][i][j];
	for (j=0;j<576;j++) xr[1][i][j]    = xr[0][i][j];
	for (j=0;j<576;j++) xr[0][i][j]    = xr_buff2[i][j];
    }
}
