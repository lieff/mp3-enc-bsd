/**********************************************************************
 * ISO MPEG Audio Subgroup Software Simulation Group (1996)
 * ISO 13818-3 MPEG-2 Audio Encoder - Lower Sampling Frequency Extension
 *
 * $Id: mdct.h,v 1.1 1996/02/14 04:04:23 rowlands Exp $
 *
 * $Log: mdct.h,v $
 * Revision 1.1  1996/02/14 04:04:23  rowlands
 * Initial revision
 *
 * Received from Mike Coleman
 **********************************************************************/

#ifndef MDCT_DOT_H
#define MDCT_DOT_H
void mdct(double *in, double *out, int block_type);
void inv_mdct(double *in, double *out, int block_type);

typedef double D32_18[SBLIMIT][18];
typedef double L3SBS[2][3][18][SBLIMIT]; /* [gr][ch] */

void mdct_sub(L3SBS (*sb_sample), double (*mdct_freq)[2][576], int stereo, III_side_info_t *l3_side, int mode_gr );
void mdct_sub_dec(double (*mdct_freq)[2][576], double inv_mdct_dec[3][2][18][32], int stereo, III_side_info_t *l3_side);
void delay(double (*xr)[2][576], int stereo);
#endif
