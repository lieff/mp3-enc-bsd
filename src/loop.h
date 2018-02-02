/**********************************************************************
 * ISO MPEG Audio Subgroup Software Simulation Group (1996)
 * ISO 13818-3 MPEG-2 Audio Encoder - Lower Sampling Frequency Extension
 *
 * $Id: loop.h,v 1.1 1996/02/14 04:04:23 rowlands Exp $
 *
 * $Log: loop.h,v $
 * Revision 1.1  1996/02/14 04:04:23  rowlands
 * Initial revision
 *
 * Received from Mike Coleman
 **********************************************************************/

#ifndef LOOP_DOT_H
#define LOOP_DOT_H
#include "common.h"
#include "l3side.h"

/**********************************************************************
 *   date   programmers                comment                        *
 * 25. 6.92  Toshiyuki Ishino          Ver 1.0                        *
 * 29.10.92  Masahiro Iwadare          Ver 2.0                        *
 * 17. 4.93  Masahiro Iwadare          Updated for IS Modification    *
 *                                                                    *
 *********************************************************************/

extern int cont_flag;

#define e              2.71828182845

/*#define SBLIMIT       32*/
#define CBLIMIT       21

#define SFB_LMAX 22
#define SFB_SMAX 13

extern int pretab[];

struct scalefac_struct
{
   int l[23];
   int s[14];
};

extern struct scalefac_struct sfBandIndex[];  /* Table B.8 -- in loop.c */


void iteration_loop( double pe[][2], double xr_org[2][2][576], III_psy_ratio *ratio,
		     III_side_info_t *l3_side, int l3_enc[2][2][576], int mean_bits,
		     int stereo, double xr_dec[2][2][576],
		     III_scalefac_t *scalefac, frame_params *fr_ps,
		     int ancillary_pad, int bitsPerFrame );


int nint( double in );

/* #define PI 3.1415926535 */
#define maximum(A,B) ( (A) > (B) ? (A) : (B) )
#define minimum(A,B) ( (A) < (B) ? (A) : (B) )
#define signum( A ) ( (A) > 0 ? 1 : -1 )

/* GLOBALE VARIABLE */

/*extern FILE     *debp;
extern FILE     *huffcp,*huffdp;
extern FILE     *musicin,*cod_music;*/

/* Beachte:Partitonen nur fuer 48 kHz */
/* andere ev. fehlerhaft  */

/* static int     scalefac_band_long[22];
  static int     scalefac_band_short[13];
*/
/*static int     huffman_tab_quad[2][16][2];
extern int     bigv_cod_tab[17][2][16][16];
extern int     bigv_dec_tab[17][512][5];
extern int     bigv_root_node[17];
extern int     count1_root_node[2];*/
/* extern int     cod_tab_info[2][31]; */




/* in max spalte war um 1 zuviel */
/* [0][i] : maximalwert der Tabelle mit index i
   [1][i] : anzahl der linbits der TAb. i        */

extern int bit_buffer[50000];

#endif
