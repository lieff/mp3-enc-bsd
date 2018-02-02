/**********************************************************************
 * ISO MPEG Audio Subgroup Software Simulation Group (1996)
 * ISO 13818-3 MPEG-2 Audio Encoder - Lower Sampling Frequency Extension
 *
 * $Id: loop.c,v 1.2 1997/01/19 22:28:29 rowlands Exp $ 
 *
 * $Log: loop.c,v $
 * Revision 1.2  1997/01/19 22:28:29  rowlands
 * Layer 3 bug fixes from Seymour Shlien
 *
 * Revision 1.1  1996/02/14 04:04:23  rowlands
 * Initial revision
 *
 * Received from Mike Coleman
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "l3side.h"
#include "loop.h"
#include "huffman.h"
#include "l3bitstream.h"
#include "reservoir.h"
#include "loop-pvt.h"
#include "sqrttab.h"

#include "pow_nint.h"
#include "huffcode.h"

/* #define DEBUG */
/* #define DEBUGSC */
/* #define PERFORM 3 */
/* If PERFORM is defined to some number, then a file encode.log
   is preduced showing the intermediate results produced by the
   outer_loop code for the frame number = PERFORM.
*/
#define BIN_SEARCH 
/* for speeding up the iteration_loop algorithm */

#ifdef PERFORM
FILE *log_output;
extern int frameNum;
float worst_xfsf_to_xmin_ratio(III_psy_xmin *l3_xmin, double xfsf[4][CBLIMIT]
                             ,int block_type,int gr,int ch);
#endif

/* New SS 20-12-96 */
#ifdef BIN_SEARCH 
int bin_search_StepSize(int desired_rate, double start, int ix[576],
           double xrs[576], gr_info * cod_info);
int count_bits();
float worst_xfsf_to_xmin_ratio();
#endif

/*
  Here are MPEG1 Table B.8 and MPEG2 Table B.1
  -- Layer III scalefactor bands.
  Index into this using a method such as:
    idx  = fr_ps->header->sampling_frequency
           + (fr_ps->header->version * 3)
*/

struct scalefac_struct sfBandIndex[6] =
{

  { /* Table B.2.b: 22.05 kHz */
    {0,6,12,18,24,30,36,44,54,66,80,96,116,140,168,200,238,284,336,396,464,522,576},
    {0,4,8,12,18,24,32,42,56,74,100,132,174,192}
  },
  { /* Table B.2.c: 24 kHz */
    {0,6,12,18,24,30,36,44,54,66,80,96,114,136,162,194,232,278,330,394,464,540,576},
    {0,4,8,12,18,26,36,48,62,80,104,136,180,192}
  },
  { /* Table B.2.a: 16 kHz */
    {0,6,12,18,24,30,36,44,45,66,80,96,116,140,168,200,238,248,336,396,464,522,576},
    {0,4,8,12,18,26,36,48,62,80,104,134,174,192}
  },
  { /* Table B.8.b: 44.1 kHz */
    {0,4,8,12,16,20,24,30,36,44,52,62,74,90,110,134,162,196,238,288,342,418,576},
    {0,4,8,12,16,22,30,40,52,66,84,106,136,192}
  },
  { /* Table B.8.c: 48 kHz */
    {0,4,8,12,16,20,24,30,36,42,50,60,72,88,106,128,156,190,230,276,330,384,576},
    {0,4,8,12,16,22,28,38,50,64,80,100,126,192}
  },
  { /* Table B.8.a: 32 kHz */
    {0,4,8,12,16,20,24,30,36,44,54,66,82,102,126,156,194,240,296,364,448,550,576},
    {0,4,8,12,16,22,30,42,58,78,104,138,180,192}
  }
};

/*
  The following table is used to implement the scalefactor
  partitioning for MPEG2 as described in section
  2.4.3.2 of the IS. The indexing corresponds to the
  way the tables are presented in the IS:

  [table_number][row_in_table][column of nr_of_sfb]
*/
static unsigned nr_of_sfb_block[6][3][4] =
{
  {
    {6, 5, 5, 5},
    {9, 9, 9, 9},
    {6, 9, 9, 9}
  },
  {
    {6, 5, 7, 3},
    {9, 9, 12, 6},
    {6, 9, 12, 6}
  },
  {
    {11, 10, 0, 0},
    {18, 18, 0, 0},
    {15,18,0,0}
  },
  {
    {7, 7, 7, 0},
    {12, 12, 12, 0},
    {6, 15, 12, 0}
  },
  {
    {6, 6, 6, 3},
    {12, 9, 9, 6},
    {6, 12, 9, 6}
  },
  {
    {8, 8, 5, 0},
    {15,12,9,0},
    {6,18,9,0}
  }
};

/*
  table of largest scalefactors for MPEG2
*/
static unsigned max_sfac_tab[6][4] =
{
    {4, 4, 3, 3},
    {4, 4, 3, 0},
    {3, 2, 0, 0},
    {4, 5, 5, 0},
    {3, 3, 3, 0},
    {2, 2, 0, 0}
};

/* Table B.6: layer3 preemphasis */
int  pretab[21] =
{
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 1, 1, 1, 2, 2, 3, 3, 3, 2
};

/* This is the scfsi_band table from 2.4.2.7 of the IS */
int scfsi_band_long[5] = { 0, 6, 11, 16, 21 };

int *scalefac_band_long  = &sfBandIndex[3].l[0];
int *scalefac_band_short = &sfBandIndex[3].s[0];



/*************************************************************************/
/*            count_bit                                                  */
/*************************************************************************/

/*
 Function: Count the number of bits necessary to code the subregion. 
*/

static inline int count_bit( int ix[576], unsigned int start, unsigned int end, unsigned int table )
{
     unsigned            linbits, ylen;
     register int        i, sum;
     register int        x,y;
     struct huffcodetab *h;
 
     if(table==0) return 0;
  
     h   = &(ht[table]);
     sum = 0;
 
     ylen    = h->ylen;
     linbits = h->linbits;
 
     if(table>15)
     { /* ESC-table is used */
         for(i=start;i<end;i+=2)
         {
             x = ix[i];
             y = ix[i+1];
             if(x>14)
             {
                 x = 15;
                 sum += linbits;
             }
             if(y>14)
             {
                 y = 15;
                 sum += linbits;
             }
 
             sum += h->hlen[(x*ylen)+y];

             if(x!=0) sum++;
             if(y!=0) sum++;
         }
     }
     else
     { /* No ESC-words */
         for(i=start;i<end;i+=2)
         {
             x = ix[i];
             y = ix[i+1];
 
             sum  += h->hlen[(x*ylen)+y];
 
             if(x!=0) sum++;
             if(y!=0) sum++;
         }
}

    return sum;
}



/************************************************************************/
/*  iteration_loop()                                                    */
/************************************************************************/
void
iteration_loop( double pe[][2], double xr_org[2][2][576], III_psy_ratio *ratio,
		III_side_info_t *l3_side, int l3_enc[2][2][576],
		int mean_bits, int stereo, double xr_dec[2][2][576],
		III_scalefac_t *scalefac, frame_params *fr_ps,
		int ancillary_pad, int bitsPerFrame ) 
{
    static int firstcall = 1;
    III_psy_xmin l3_xmin;
    gr_info *cod_info;
    layer *info;
    int *main_data_begin;

    int max_bits;
    int ch, gr, sfb, i, mode_gr;
    static int huffman_read_flag = 0; 

    double xr[2][2][576];
    I576  *ix;

    main_data_begin = &l3_side->main_data_begin;
    l3_side->resvDrain = 0;

    if ( firstcall )
    {
	*main_data_begin = 0;
	firstcall = 0;
#ifdef PERFORM
	log_output = fopen("encode.log","w");
#endif
    }

    info = fr_ps->header;
    mode_gr = (info->version == 1) ? 2 : 1;

    scalefac_band_long  = &sfBandIndex[info->sampling_frequency + (info->version * 3)].l[0];
    scalefac_band_short = &sfBandIndex[info->sampling_frequency + (info->version * 3)].s[0];

    /* reading huffman code table */
#if 0 /* tables now incorporated in exe */
    if (huffman_read_flag == 0) {
        FILE *fi = OpenTableFile( "huffcode" );
        if ( fi == NULL )
            exit( EXIT_FAILURE );
        read_huffcodetab( fi );
        huffman_read_flag++;
        fclose( fi );
    }
#endif

    for ( gr = 0; gr < mode_gr; gr++ )
    {
        for ( ch = 0; ch < stereo; ch++ )
	{
            for ( i = 0; i < 576; i++ ) 
                xr[gr][ch][i] = xr_org[gr][ch][i];
	}
    }

    ResvFrameBegin( fr_ps, l3_side, mean_bits, bitsPerFrame );

    for ( gr = 0; gr < mode_gr; gr++ )
    {
        for ( ch = 0; ch < stereo; ch++ )
        {
            ix = (I576 *) l3_enc[gr][ch];
            cod_info = (gr_info *) &(l3_side->gr[gr].ch[ch]);
            gr_deco(cod_info);
            calc_xmin( xr, ratio, cod_info, &l3_xmin, gr, ch );
	    
#ifdef DEBUG
            printf( "----- start gr[%d] ch[%1d] : block_type=%1d, window_switching_flag=%1d (loop)\n",
                    gr,ch, cod_info->block_type, cod_info->window_switching_flag ); 
#endif
	    
            if ( info->version == 1 )
                calc_scfsi( xr[gr][ch], l3_side, &l3_xmin, ch, gr );
	    
            /* calculation of number of available bit( per granule ) */
	    max_bits = ResvMaxBits( fr_ps, l3_side, &pe[gr][ch], mean_bits );
#ifdef DEBUG
            printf( " max_bits = %d, mean_bits = %d (iteration_loop)\n", max_bits, mean_bits ); 
#endif
	    
            /* reset of iteration variables */
	    
            for ( sfb = 0; sfb < 21; sfb++ )
                scalefac->l[gr][ch][sfb] = 0;
            for ( sfb = 0; sfb < 13; sfb++ )
                for ( i = 0; i < 3; i++ )
                    scalefac->s[gr][ch][sfb][i] = 0;

	    for ( i = 0; i < 4; i++ )
		cod_info->slen[i] = 0;
	    cod_info->sfb_partition_table = &nr_of_sfb_block[0][0][0];

            cod_info->part2_3_length    = 0;
            cod_info->big_values        = 0;
            cod_info->count1            = 0;
            cod_info->scalefac_compress = 0;
            cod_info->table_select[0]   = 0;
            cod_info->table_select[1]   = 0;
            cod_info->table_select[2]   = 0;
            cod_info->subblock_gain[0]  = 0;
            cod_info->subblock_gain[1]  = 0;
            cod_info->subblock_gain[2]  = 0;
            cod_info->region0_count     = 0;
            cod_info->region1_count     = 0;
            cod_info->part2_length      = 0;
            cod_info->preflag           = 0;
            cod_info->scalefac_scale    = 0;
            cod_info->quantizerStepSize = 0.0;
            cod_info->count1table_select= 0;
            
            /* all spectral values zero ? */
            if ( fabs(xr_max(xr[gr][ch], 0, 576)) != 0.0 )
            {
                cod_info->quantizerStepSize =
			 (double) quantanf_init( xr[gr][ch] );
                cod_info->part2_3_length = outer_loop( xr, max_bits, &l3_xmin,
                                                       l3_enc, fr_ps, scalefac,
                                                       gr, ch, l3_side );
            }
	    ResvAdjust( fr_ps, cod_info, l3_side, mean_bits );

	    cod_info->global_gain = nint( cod_info->quantizerStepSize + 210.0 );
	    assert( cod_info->global_gain < 256 );
        } /* for ch */
    } /* for gr */
    ResvFrameEnd( fr_ps, l3_side, mean_bits );
}



/************************************************************************/
/*  quantanf_init                                                       */
/************************************************************************/
int quantanf_init( double xr[576] )
/* Function: Calculate the first quantization step quantanf.       */
{
    int i, tp = 0;
    double system_const, minlimit;
    double sfm = 0.0, sum1 = 0.0, sum2 = 0.0;
    
    system_const = 8.0;
    minlimit = -100.0;

    for ( i = 0; i < 576; i++ )
    {
        if ( xr[i] != 0 )
	{
            double tpd = xr[i] * xr[i];
            sum1 += log( tpd );
            sum2 += tpd;
        }
    }
    if ( sum2 != 0.0 )
    {
        sfm = exp( sum1 / 576.0 ) / (sum2 / 576.0);
        tp = nint( system_const * log(sfm) );
	if ( tp < minlimit )
	    tp = minlimit;
#ifdef DEBUG
        printf(" quantanf = %d (quantanf_init)\n",tp );
#endif
    }
      return(tp-70.0); /* SS 19-12-96. Starting value of
                          global_gain or quantizerStepSize 
                          has to be reduced for iteration_loop
                       */
}





/************************************************************************/
/*  outer_loop                                                          */
/************************************************************************/
/*  Function: The outer iteration loop controls the masking conditions  */
/*  of all scalefactorbands. It computes the best scalefac and          */
/*  global gain. This module calls the inner iteration loop             */
/************************************************************************/
int outer_loop(
    double xr[2][2][576],     /*  magnitudes of the spectral values */
    int max_bits,
    III_psy_xmin  *l3_xmin,   /* the allowed distortion of the scalefactor */
    int l3_enc[2][2][576],    /* vector of quantized values ix(0..575) */
    frame_params *fr_ps,
    III_scalefac_t *scalefac, /* scalefactors */
    int gr, int ch, III_side_info_t *l3_side )
{
    int status ;
    int scalesave_l[CBLIMIT], scalesave_s[CBLIMIT][3];
    int sfb, bits, huff_bits, save_preflag, save_compress;
    double xfsf[4][CBLIMIT];
    int i, over, iteration;
    float max_ratio;


/* D576 *xrs; */ /* to eliminate warning messages from gcc compiler */
/* I576 *ix; */  /* replace this code with below. S. Shlien 15-1-97 */

  double *xrs; 
  int *ix;  
  gr_info *cod_info = &l3_side->gr[gr].ch[ch].tt;

/* xrs = (D576 *) &xr[gr][ch][0]; */ 
/* ix  = (I576 *) l3_enc[gr][ch]; */
 

xrs = (double *) &(xr[gr][ch][0]); 
ix  = (int *) &(l3_enc[gr][ch][0]);


    iteration = 0;
#ifdef PERFORM
    if(frameNum == PERFORM)
    fprintf(log_output,"\n\nframe = %d ch = %d gr= %d\n",frameNum,ch,gr);
#endif
    do 
    {
	iteration += 1;
	cod_info->part2_length = part2_length( scalefac, fr_ps, gr, ch, l3_side );
        huff_bits = max_bits - cod_info->part2_length;

#ifdef BIN_SEARCH
	if(iteration == 1)
        {
       bin_search_StepSize(max_bits,cod_info->quantizerStepSize,
          ix,xrs,cod_info); /* speeds things up a bit */
        }
#endif
#ifdef PERFORM
	if(frameNum==PERFORM)
          fprintf(log_output,"\n    Interim Results %d\n\n",iteration);
        bits = test_inner_loop( xr, l3_enc, huff_bits, cod_info, gr, ch,
         xfsf,l3_xmin);
#else
        bits = inner_loop( xr, l3_enc, huff_bits, cod_info, gr, ch );
#endif

        calc_noise( &xr[gr][ch][0], &l3_enc[gr][ch][0], cod_info, xfsf ); /* distortion calculation */

        for ( sfb = 0; sfb < CBLIMIT; sfb++ ) /* save scaling factors */
            scalesave_l[sfb] = scalefac->l[gr][ch][sfb];

        for ( sfb = 0; sfb < SFB_SMAX; sfb++ )
            for ( i = 0; i < 3; i++ )
                scalesave_s[sfb][i] = scalefac->s[gr][ch][sfb][i];
        
        save_preflag  = cod_info->preflag;
        save_compress = cod_info->scalefac_compress;

        preemphasis( &xr[gr][ch][0], xfsf, l3_xmin, gr, ch, l3_side );



#ifdef PERFORM 
    if(frameNum == PERFORM)  
    {
    fprintf(log_output,"\nbits = %d  huff_bits= %d ", bits,huff_bits);
    fprintf(log_output," max_bits = %d\n",max_bits);
    fprintf(log_output,"Stepsize = %f ",cod_info->quantizerStepSize);
    fprintf(log_output," scale_bits    = %d \n", cod_info->part2_length );  
    print_scalefacs(log_output,scalefac,cod_info->block_type,gr,ch);
    /*if (gr==0 && ch==0)
      print_quantized_values(log_output,l3_enc[gr][ch] ,cod_info);*/
    max_ratio = worst_xfsf_to_xmin_ratio(l3_xmin,xfsf,cod_info->block_type,gr,ch);
    fprintf(log_output,"max_ratio = %6.2f\n",max_ratio );
    print_ratios(log_output,l3_xmin,xfsf,cod_info->block_type,gr,ch);
    fprintf(log_output,"\n\n");
    fflush(log_output);
    }
#endif
        over = amp_scalefac_bands( &xr[gr][ch][0], xfsf, l3_xmin,
                                   l3_side, scalefac, gr, ch, iteration );

#if 1
        if ( (status = loop_break(scalefac, cod_info, gr, ch)) == 0 )
	{
	    if ( fr_ps->header->version == 1 )
		status = scale_bitcount( scalefac, cod_info, gr, ch );
	    else
		status = scale_bitcount_lsf( scalefac, cod_info, gr, ch );
	}
#else
        status = loop_break( scalefac, cod_info, gr, ch );
	if ( fr_ps->header->version == 1 )
	    status += scale_bitcount( scalefac, cod_info, gr, ch );
	else
	    status += scale_bitcount_lsf( scalefac, cod_info, gr, ch );
#endif

    }
    while ( (status == 0) && (over > 0) );

    cod_info->preflag = save_preflag;
    cod_info->scalefac_compress = save_compress;

    for ( sfb = 0; sfb < 21; sfb++ )
        scalefac->l[gr][ch][sfb] = scalesave_l[sfb];    

    for ( i = 0; i < 3; i++ )
        for ( sfb = 0; sfb < 12; sfb++ )
            scalefac->s[gr][ch][sfb][i] = scalesave_s[sfb][i];    

    cod_info->part2_length   = part2_length( scalefac, fr_ps, gr, ch, l3_side );
    cod_info->part2_3_length = cod_info->part2_length + bits;

#ifdef PERFORM 
    if(frameNum == PERFORM)
    {
    fprintf(log_output,"\n  Final Results\n");
    fprintf(log_output,"bits = %d  huff_bits= %d", bits,huff_bits);
    fprintf(log_output," max_bits = %d\n",max_bits);
    fprintf(log_output," Stepsize = %f ",cod_info->quantizerStepSize);
    fprintf(log_output, " scale_bits    = %d \n", cod_info->part2_length );  
    max_ratio = worst_xfsf_to_xmin_ratio(l3_xmin,xfsf,cod_info->block_type,gr,ch);
    print_scalefacs(log_output,scalefac,cod_info->block_type,gr,ch);
    fprintf(log_output,"max_ratio = %6.2f\n",max_ratio );
    print_ratios(log_output,l3_xmin,xfsf,cod_info->block_type,gr,ch);
    fflush(log_output);
    }
#endif    
    return cod_info->part2_3_length;
}




/***************************************************************************/ 
/*         inner_loop                                                      */ 
/***************************************************************************/ 
/* The code selects the best quantizerStepSize for a particular set
/* of scalefacs                                                            */
 
int
inner_loop( double xr[2][2][576], int l3_enc[2][2][576], int max_bits,
	    gr_info *cod_info, int gr, int ch )
{
    int bits, c1bits, bvbits;
    double *xrs;  /*  D576  *xr; */
    int *ix;  /*  I576  *ix; */
    xrs = &xr[gr][ch][0];
    ix = l3_enc[gr][ch];

    assert( max_bits >= 0 );
    cod_info->quantizerStepSize -= 1.0;;
    do
    {
        do
        {
            cod_info->quantizerStepSize += 1.0;
            quantize( xrs, ix, cod_info );
        }
        while ( ix_max(ix, 0, 576) > 8191 + 14 ); /* within table range? */

        calc_runlen( ix, cod_info );  /*rzero,count1,big_values*/
        bits = c1bits = count1_bitcount( ix, cod_info );  /*count1_table selection*/
        subdivide( cod_info );  /* bigvalues sfb division */
        bigv_tab_select( ix, cod_info );  /* codebook selection*/
        bits += bvbits = bigv_bitcount( ix, cod_info );  /* bit count */

#ifdef PERFORM
        if(frameNum == PERFORM)
	fprintf(log_output,"StepSize=%f bits = %d huff_bits = %d\n",
           cod_info->quantizerStepSize,bits,max_bits);
#endif

    }
    while ( bits > max_bits );

    return bits;
}



/***************************************************************************/ 
/*        calc_scfsi                                                       */ 
/***************************************************************************/ 
/* calculation of the scalefactor select information ( scfsi )        */

void calc_scfsi( double  xr[576], III_side_info_t *l3_side,
	    III_psy_xmin *l3_xmin, int ch, int gr )
{
    static int en_tot[2][2]; /* ch,gr */
    static int en[2][2][21];
    static int xm[2][2][21];
    static int xrmax[2][2];

    int en_tot_krit        = 10;
    int en_dif_krit        = 100;
    int en_scfsi_band_krit = 10;
    int xm_scfsi_band_krit = 10;

    int scfsi_band;
    unsigned scfsi_set;

    int sfb, start, end, i;
    int condition = 0;
    double temp, log2 = log( 2.0 );
    gr_info *cod_info = &l3_side->gr[gr].ch[ch].tt;

    xrmax[gr][ch] = xr_max( xr, 0, 576 );
    scfsi_set = 0;

    /* the total energy of the granule */    
    for ( temp = 0.0, i = 0; i < 576; i++ )
        temp += xr[i] * xr[i];
    if ( temp == 0.0 )
        en_tot[gr][ch] = 0.0;
    else
        en_tot[gr][ch] = log( temp ) / log2 ;

    /* the energy of each scalefactor band, en */
    /* the allowed distortion of each scalefactor band, xm */

    if ( cod_info->window_switching_flag == 0 ||
         cod_info->block_type != 2 )
        for ( sfb = 0; sfb < 21; sfb++ )
        {
            start = scalefac_band_long[ sfb ];
            end   = scalefac_band_long[ sfb+1 ];

            for ( temp = 0.0, i = start; i < end; i++ )
                temp += xr[i] * xr[i];
            if ( temp == 0.0 )
                en[gr][ch][sfb] = 0.0;
            else
                en[gr][ch][sfb] = log( temp )/ log2;

            if ( l3_xmin->l[gr][ch][sfb] == 0.0 )
                xm[gr][ch][sfb] = 0.0;
            else
                xm[gr][ch][sfb] = log( l3_xmin->l[gr][ch][sfb] ) / log2;
        }
    if ( gr == 1 )
    {
        int gr2, tp;

        for ( gr2 = 0; gr2 < 2; gr2++ )
        {
            /* The spectral values are not all zero */
            if ( xrmax[ch][gr2] != 0.0 )
                condition++;
            /* None of the granules contains short blocks */
            if ( (cod_info->window_switching_flag == 0) ||
                 (cod_info->block_type != 2) )
                condition++;
        }
        if ( abs(en_tot[0] - en_tot[1]) < en_tot_krit )
            condition++;
        for ( tp = 0, sfb = 0; sfb < 21; sfb++ ) 
            tp += abs( en[ch][0][sfb] - en[ch][1][sfb] );
        if ( tp < en_dif_krit ) 
            condition++;

        if ( condition == 6 )
        {
            for ( scfsi_band = 0; scfsi_band < 4; scfsi_band++ )
            {
                int sum0 = 0, sum1 = 0;
                l3_side->scfsi[ch][scfsi_band] = 0;
                start = scfsi_band_long[scfsi_band];
                end   = scfsi_band_long[scfsi_band+1];
                for ( sfb = start; sfb < end; sfb++ )
                { 
                    sum0 += abs( en[ch][0][sfb] - en[ch][1][sfb] );
                    sum1 += abs( xm[ch][0][sfb] - xm[ch][1][sfb] );
                }

                if ( sum0 < en_scfsi_band_krit && sum1 < xm_scfsi_band_krit )
		{
                    l3_side->scfsi[ch][scfsi_band] = 1;
		    scfsi_set |= (1 << scfsi_band);
		}
                else
                    l3_side->scfsi[ch][scfsi_band] = 0;
            } /* for scfsi_band */
#ifdef DEBUG
	    fprintf( stderr, "calc_scfsi: scfsi_set = 0x%02x\n", scfsi_set );
#endif
        } /* if condition == 6 */
        else
            for ( scfsi_band = 0; scfsi_band < 4; scfsi_band++ )
                l3_side->scfsi[ch][scfsi_band] = 0;
    } /* if gr == 1 */
}



/***************************************************************************/ 
/*        part2_length                                                     */ 
/***************************************************************************/ 

/* calculates the number of bits needed to encode the scalefacs in the     */
/* main data block                                                         */

int part2_length( III_scalefac_t *scalefac, frame_params *fr_ps,
	      int gr, int ch, III_side_info_t *si )
{
    int slen1, slen2, slen3, slen4, bits, sfb, window, partition;
    gr_info *gi = &si->gr[gr].ch[ch].tt;

    bits = 0;
    if ( fr_ps->header->version == 1 )
    {
	static int slen1_tab[16] = { 0, 0, 0, 0, 3, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4 };
	static int slen2_tab[16] = { 0, 1, 2, 3, 0, 1, 2, 3, 1, 2, 3, 1, 2, 3, 2, 3 };

	slen1 = slen1_tab[ gi->scalefac_compress ];
	slen2 = slen2_tab[ gi->scalefac_compress ];

	if ( (gi->window_switching_flag == 1) && (gi->block_type == 2) )
	{
	    if ( gi->mixed_block_flag )
	    {
		bits += (8 * slen1) + (9 * slen1) + (18 * slen2);
	    }
	    else
	    {
		bits += (18 * slen1) + (18 * slen2);
	    }
	}
	else
	{
	    if ( (gr == 0) || (si->scfsi[ch][0] == 0) )
		bits += (6 * slen1);

	    if ( (gr == 0) || (si->scfsi[ch][1] == 0) )
		/*  bits += (6 * slen1);  This is wrong SS 19-12-96 */
		bits += (5 * slen1);

	    if ( (gr == 0) || (si->scfsi[ch][2] == 0) )
		/*  bits += (6 * slen1);   This is wrong SS 19-12-96 */
		bits += (5 * slen2);

	    if ( (gr == 0) || (si->scfsi[ch][3] == 0) )
		/* bits += (6 * slen1);   This is wrong SS 19-12-96 */
		bits += (5 * slen2);
	}
    }
    else
    {   /* MPEG 2 */
	assert( gi->sfb_partition_table );
	for ( partition = 0; partition < 4; partition++ )
	    bits += gi->slen[partition] * gi->sfb_partition_table[partition];
    }
    return bits;
}



/*************************************************************************/
/*            scale_bitcount                                             */
/*************************************************************************/

/* Also calculates the number of bits necessary to code the scalefactors. */

int scale_bitcount( III_scalefac_t *scalefac, gr_info *cod_info,
		int gr, int ch )
{
    int i, k, sfb, max_slen1 = 0, max_slen2 = 0, /*a, b, */ ep = 2;

    static int slen1[16] = { 0, 0, 0, 0, 3, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4 };
    static int slen2[16] = { 0, 1, 2, 3, 0, 1, 2, 3, 1, 2, 3, 1, 2, 3, 2, 3 };
    static int pow2[5]   = { 1, 2, 4, 8, 16 };

    if ( cod_info->window_switching_flag != 0 && cod_info->block_type == 2 )
    {
        if ( cod_info->mixed_block_flag == 0 ) 
        {
            /* a = 18; b = 18;  */
            for ( i = 0; i < 3; i++ )
            {
                for ( sfb = 0; sfb < 6; sfb++ )
                    if ( scalefac->s[gr][ch][sfb][i] > max_slen1 )
                        max_slen1 = scalefac->s[gr][ch][sfb][i];
                for (sfb = 6; sfb < 12; sfb++ )
                    if ( scalefac->s[gr][ch][sfb][i] > max_slen2 )
                        max_slen2 = scalefac->s[gr][ch][sfb][i];
            }
        }
        else
        {/* mixed_block_flag = 1 */
            /* a = 17; b = 18;  */
            for ( sfb = 0; sfb < 8; sfb++ )
                if ( scalefac->l[gr][ch][sfb] > max_slen1 )
                    max_slen1 = scalefac->l[gr][ch][sfb];
            for ( i = 0; i < 3; i++ )
            {
                for ( sfb = 3; sfb < 6; sfb++ )
                    if ( scalefac->s[gr][ch][sfb][i] > max_slen1 )
                        max_slen1 = scalefac->s[gr][ch][sfb][i];
                for ( sfb = 6; sfb < 12; sfb++ )
                    if ( scalefac->s[gr][ch][sfb][i] > max_slen2 )
                        max_slen2 = scalefac->s[gr][ch][sfb][i];
            }
        }
    }
    else
    { /* block_type == 1,2,or 3 */
        /* a = 11; b = 10;   */
        for ( sfb = 0; sfb < 11; sfb++ )
            if ( scalefac->l[gr][ch][sfb] > max_slen1 )
                max_slen1 = scalefac->l[gr][ch][sfb];
        for ( sfb = 11; sfb < 21; sfb++ )
            if ( scalefac->l[gr][ch][sfb] > max_slen2 )
                max_slen2 = scalefac->l[gr][ch][sfb];
    }

    for ( k = 0; k < 16; k++ )
    {
        if ( (max_slen1 < pow2[slen1[k]]) && (max_slen2 < pow2[slen2[k]]) )
        { 
            ep = 0;
            break;
        } 
    }

    if ( ep == 0 )
        cod_info->scalefac_compress = k;
    return ep;
}




/*************************************************************************/
/*            scale_bitcount_lsf                                         */
/*************************************************************************/

/* Also counts the number of bits to encode the scalefacs but for MPEG 2 */ 
/* Lower sampling frequencies  (24, 22.05 and 16 kHz.)                   */
 
/*  This is reverse-engineered from section 2.4.3.2 of the MPEG2 IS,     */
/* "Audio Decoding Layer III"                                            */

int scale_bitcount_lsf( III_scalefac_t *scalefac, gr_info *cod_info,
		    int gr, int ch )
{
    int table_number, row_in_table, partition, nr_sfb, window, over;
    int i, k, sfb, max_sfac[ 4 ];
    unsigned *partition_table;

    /*
      Set partition table. Note that should try to use table one,
      but do not yet...
    */
    if ( cod_info->preflag )
	table_number = 2;
    else
	table_number = 0;

    for ( i = 0; i < 4; i++ )
	max_sfac[i] = 0;

    if ( cod_info->window_switching_flag != 0 && cod_info->block_type == 2 )
    {
        if ( cod_info->mixed_block_flag == 0 ) 
        {
	    row_in_table = 1;
	    partition_table = &nr_of_sfb_block[table_number][row_in_table][0];
	    for ( sfb = 0, partition = 0; partition < 4; partition++ )
	    {
		nr_sfb = partition_table[ partition ] / 3;
		for ( i = 0; i < nr_sfb; i++, sfb++ )
		    for ( window = 0; window < 3; window++ )
			if ( scalefac->s[gr][ch][sfb][window] > max_sfac[partition] )
			    max_sfac[partition] = scalefac->s[gr][ch][sfb][window];
	    }
        }
        else
        {/* mixed_block_flag = 1 */
	    row_in_table = 2;
	    partition_table = &nr_of_sfb_block[table_number][row_in_table][0];
	    partition = 0;
	    nr_sfb = partition_table[ partition ];
	    for ( i = 0; i < nr_sfb; i++, sfb++ )
		if ( scalefac->l[gr][ch][sfb] > max_sfac[partition] )
		    max_sfac[partition] = scalefac->l[gr][ch][sfb];
	    
	    for ( sfb = 0, partition = 1; partition < 4; partition++ )
	    {
		nr_sfb = partition_table[ partition ] / 3;
		for ( i = 0; i < nr_sfb; i++, sfb++ )
		    for ( window = 0; window < 3; window++ )
			if ( scalefac->s[gr][ch][sfb][window] > max_sfac[partition] )
			    max_sfac[partition] = scalefac->s[gr][ch][sfb][window];
	    }
        }
    }
    else
    {
	row_in_table = 0;
	partition_table = &nr_of_sfb_block[table_number][row_in_table][0];
	partition = 0;
	for ( sfb = 0, partition = 0; partition < 4; partition++ )
	{
	    nr_sfb = partition_table[ partition ];
	    for ( i = 0; i < nr_sfb; i++, sfb++ )
		if ( scalefac->l[gr][ch][sfb] > max_sfac[partition] )
		    max_sfac[partition] = scalefac->l[gr][ch][sfb];
	}
    }

    for ( over = 0, partition = 0; partition < 4; partition++ )
    {
	if ( max_sfac[partition] > max_sfac_tab[table_number][partition] )
	    over++;
    }
    if ( !over )
    {
	/*
	  Since no bands have been over-amplified, we can set scalefac_compress
	  and slen[] for the formatter
	*/
	static int log2tab[] = { 0, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4 };

	unsigned slen1, slen2, slen3, slen4;

        cod_info->sfb_partition_table = &nr_of_sfb_block[table_number][row_in_table][0];
	for ( partition = 0; partition < 4; partition++ )
	    cod_info->slen[partition] = log2tab[max_sfac[partition]];

	/* set scalefac_compress */
	slen1 = cod_info->slen[ 0 ];
	slen2 = cod_info->slen[ 1 ];
	slen3 = cod_info->slen[ 2 ];
	slen4 = cod_info->slen[ 3 ];

	switch ( table_number )
	{
	  case 0:
	    cod_info->scalefac_compress = (((slen1 * 5) + slen2) << 4)
		+ (slen3 << 2)
		+ slen4;
	    break;

	  case 1:
	    cod_info->scalefac_compress = 400
		+ (((slen1 * 5) + slen2) << 2)
		+ slen3;
	    break;

	  case 2:
	    cod_info->scalefac_compress = 500 + (slen1 * 3) + slen2;
	    break;

	  default:
	    fprintf( stderr, "intensity stereo not implemented yet\n" );
	    exit( EXIT_FAILURE );
	    break;
	}
    }
#ifdef DEBUG
    if ( over ) 
        printf( "---WARNING !! Amplification of some bands over limits\n" );
#endif
    return over;
}


/*************************************************************************/
/*            calc_noise                                                 */
/*************************************************************************/

/*   Function: Calculate the distortion introduced by the quantization   */
/*   in each scale factor band.                                          */

#define POW_4_3_SIZE 1024

double pow_4_3_tab[POW_4_3_SIZE];

void calc_noise( double xr[576], int ix[576], gr_info *cod_info,
	    double xfsf[4][CBLIMIT] )
{
    int start, end, sfb, l, i, ixst;
    double sum,step,bw;
    static int inited = 0;
    
    D192_3 *xr_s;
    I192_3 *ix_s;

    if (!inited) {
      inited = 1;
      for (i = 0; i < POW_4_3_SIZE; i++)
	pow_4_3_tab[i] = pow ((double) i, 4.0/3.0);
    }
    
    xr_s = (D192_3 *) xr;
    ix_s = (I192_3 *) ix;

    step = pow( 2.0, (cod_info->quantizerStepSize) * 0.25 );
    for ( sfb = 0; sfb < cod_info->sfb_lmax; sfb++ )
    {
        start = scalefac_band_long[ sfb ];
        end   = scalefac_band_long[ sfb+1 ];
	bw = end - start;

        for ( sum = 0.0, l = start; l < end; l++ )
        {
            double temp;
	    
	    ixst = ix[l];
	    if (ixst >= POW_4_3_SIZE)
	      temp = fabs( xr[l] ) - pow( (double) ixst, 4.0 / 3.0 ) * step;
	    else
	      temp = fabs( xr[l] ) - pow_4_3_tab[ixst] * step;
            sum += temp * temp; 
        }
        xfsf[0][sfb] = sum / bw;
    }

    for ( i = 0; i < 3; i++ )
    {
        step = pow( 2.0, (cod_info->quantizerStepSize) * 0.25 ); /* subblock_gain ? */
        for ( sfb = cod_info->sfb_smax; sfb < 12; sfb++ )
        {
            start = scalefac_band_short[ sfb ];
            end   = scalefac_band_short[ sfb+1 ];
	    bw = end - start;
            
            for ( sum = 0.0, l = start; l < end; l++ )
            {
                double temp;
		
		ixst = (*ix_s)[l][i];
		if (ixst >= POW_4_3_SIZE)
		  temp = fabs( (*xr_s)[l][i] ) - pow( (double) ixst, 4.0 / 3.0 ) * step;
		else
		  temp = fabs( (*xr_s)[l][i] ) - pow_4_3_tab[ixst] * step;
                sum += temp * temp;
            }       
            xfsf[i+1][sfb] = sum / bw;
        }
    }
}




/*************************************************************************/
/*            calc_xmin                                                  */
/*************************************************************************/

/*
  Calculate the allowed distortion for each scalefactor band,
  as determined by the psychoacoustic model.
  xmin(sb) = ratio(sb) * en(sb) / bw(sb)
*/

void calc_xmin( double xr[2][2][576], III_psy_ratio *ratio,
	   gr_info *cod_info, III_psy_xmin *l3_xmin,
	   int gr, int ch )
{
    int start, end, sfb, l, b;
    double en, bw;

    D192_3 *xr_s;

    xr_s = (D192_3 *) xr[gr][ch] ;

    for ( sfb = cod_info->sfb_smax; sfb < SFB_SMAX - 1; sfb++ )
    {
        start = scalefac_band_short[ sfb ];
        end   = scalefac_band_short[ sfb + 1 ];
	bw = end - start;
        for ( b = 0; b < 3; b++ )
        {
            for ( en = 0.0, l = start; l < end; l++ )
                en += (*xr_s)[l][b] * (*xr_s)[l][b];
            l3_xmin->s[gr][ch][sfb][b] = ratio->s[gr][ch][sfb][b] * en / bw;
        }
    }

    for ( sfb = 0; sfb < cod_info->sfb_lmax; sfb++ )
    {
        start = scalefac_band_long[ sfb ];
        end   = scalefac_band_long[ sfb+1 ];
	bw = end - start;

        for ( en = 0.0, l = start; l < end; l++ )
            en += xr[gr][ch][l] * xr[gr][ch][l];
        l3_xmin->l[gr][ch][sfb] = ratio->l[gr][ch][sfb] * en / bw;
    }
}



/*************************************************************************/
/*            loop_break                                                 */
/*************************************************************************/

/*  Function: Returns zero if there is a scalefac which has not been
    amplified. Otherwise it returns one. 
*/

int loop_break( III_scalefac_t *scalefac, gr_info *cod_info,
	    int gr, int ch )
{
    int i, sfb, temp = 1;

    for ( sfb = 0; sfb < cod_info->sfb_lmax; sfb++ )
        if ( scalefac->l[gr][ch][sfb] == 0 )
            temp = 0;

    for ( sfb = cod_info->sfb_smax; sfb < 12; sfb++ )
        for ( i = 0; i < 3; i++ )
            if ( scalefac->s[gr][ch][sfb][i] == 0 )
                temp = 0;
#ifdef DEBUG
    if ( temp != 0 )
        printf( "---WARNING !! All scalefactor bands amplified\n" );
#endif
    return temp;
}



/*************************************************************************/
/*            preemphasis                                                */
/*************************************************************************/

/*
  See ISO 11172-3  section  C.1.5.4.3.4
*/

void preemphasis( double xr[576], double xfsf[4][CBLIMIT],
	     III_psy_xmin  *l3_xmin,
	     int gr, int ch, III_side_info_t *l3_side )
{
    int i, sfb, start, end, scfsi_band, over;
    double ifqstep;
    gr_info *cod_info = &l3_side->gr[gr].ch[ch].tt;

    if ( gr == 1 )
    {
	/*
	  If the second granule is being coded and scfsi is active in
	  at least one scfsi_band, the preemphasis in the second granule
	  is set equal to the setting in the first granule
	*/
	for ( scfsi_band = 0; scfsi_band < 4; scfsi_band++ )
	    if ( l3_side->scfsi[ch][scfsi_band] )
	    {
		cod_info->preflag = l3_side->gr[0].ch[ch].tt.preflag;
		return;
	    }
	
    }

    /*
      Preemphasis is switched on if in all the upper four scalefactor
      bands the actual distortion exceeds the threshold after the
      first call of the inner loop
    */
    if ( cod_info->block_type != 2 && cod_info->preflag == 0 )
    {	
	over = 0;
	for ( sfb = 17; sfb < 21; sfb++ )
	    if ( xfsf[0][sfb] > l3_xmin->l[gr][ch][sfb] )
		over++;

	if (over == 4 )
	{
	    cod_info->preflag = 1;
	    ifqstep = ( cod_info->scalefac_scale == 0 ) ? sqrt(2.)
		: pow( 2.0, (0.5 * (1.0 + (double) cod_info->scalefac_scale)) );

	    for ( sfb = 0; sfb < cod_info->sfb_lmax; sfb++ )
	    {
		l3_xmin->l[gr][ch][sfb] *= pow( ifqstep, 2.0 * (double) pretab[sfb] );
		start = scalefac_band_long[ sfb ];
		end   = scalefac_band_long[ sfb+1 ];
		for( i = start; i < end; i++ )
		    xr[i] *= pow( ifqstep, (double) pretab[sfb] );
	    }
	}
    }
}


/*************************************************************************/
/*            amp_scalefac_bands                                         */
/*************************************************************************/

/* 
  Amplify the scalefactor bands that violate the masking threshold.
  See ISO 11172-3 Section C.1.5.4.3.5
*/

int amp_scalefac_bands( double xr[576], double xfsf[4][CBLIMIT],
		    III_psy_xmin *l3_xmin, III_side_info_t *l3_side,
		    III_scalefac_t *scalefac,
		    int gr, int ch, int iteration )
{
    int start, end, l, sfb, i, scfsi_band, over = 0;
    double ifqstep, ifqstep2;
    D192_3 *xr_s;
    gr_info *cod_info, *gr0;
    int copySF, preventSF;
#define PRECALC_SIZE 1024 /* WAS 256 !!! */
    static double pow43[PRECALC_SIZE];
    static int init=0;
    
    cod_info = &l3_side->gr[gr].ch[ch].tt;
    gr0      = &l3_side->gr[0].ch[ch].tt;

    if (init==0) {
      init++;
      for(i=0;i<PRECALC_SIZE;i++)
        pow43[i] = pow((double)i, 4.0/3.0);
    }
      
    xr_s = (D192_3 *) xr;
    copySF = 0;
    preventSF = 0;

    if ( cod_info->scalefac_scale == 0 )
	ifqstep = sqrt( 2.0 );
    else
	ifqstep = pow( 2.0, 0.5 * (1.0 + (double) cod_info->scalefac_scale) );

    if ( gr == 1 )
    {
	/*
	  If the second granule is being coded and scfsi is active in at
	  least one scfsi_band...
	*/
	for ( scfsi_band = 0; scfsi_band < 4; scfsi_band++ )
	    if ( l3_side->scfsi[ch][scfsi_band] )
	    {
		/*
		  a) ifqstep has to be set similar to the
		   first granule...
		*/
		if ( gr0->scalefac_scale == 0 )
		    ifqstep = sqrt( 2.0 );
		else
		    ifqstep = pow( 2.0, 0.5 * (1.0 + (double) gr0->scalefac_scale) );

		if ( iteration == 1 )
		{
		    /*
		      b) If it is the first iteration, the scalefactors
		      of scalefactor bands in which scfsi is enabled
		      must be taken from the first granule
		    */  
		    copySF = 1;
		}
		else
		{
		    /*
		      c) If it is not the first iteration, the amplification
		      must be prevented for scalefactor bands in which
		      scfsi is enabled
		    */
		    preventSF = 1;
		}
		break;
	    }
	
    }

    ifqstep2 = ifqstep * ifqstep;
    scfsi_band = 0;
    
    for ( sfb = 0; sfb < cod_info->sfb_lmax; sfb++ )
    {
	if ( copySF || preventSF )
	{
	    if ( sfb == scfsi_band_long[scfsi_band + 1] )
		scfsi_band += 1;
	    if ( l3_side->scfsi[ch][scfsi_band] )
	    {
		if ( copySF )
		    scalefac->l[gr][ch][sfb] = scalefac->l[0][ch][sfb];
		continue;
	    }
	}	    
	if ( xfsf[0][sfb] > l3_xmin->l[gr][ch][sfb] )
	{
	    over++;
	    l3_xmin->l[gr][ch][sfb] *= ifqstep2;
	    scalefac->l[gr][ch][sfb]++;
	    start = scalefac_band_long[sfb];
	    end   = scalefac_band_long[sfb+1];
	    for ( l = start; l < end; l++ )
		xr[l] *= ifqstep;
	}
    }

    /*
      Note that scfsi is not enabled for frames containing
      short blocks
    */
    for ( i = 0; i < 3; i++ )
        for ( sfb = cod_info->sfb_smax; sfb < 12; sfb++ )
            if ( xfsf[i+1][sfb] > l3_xmin->s[gr][ch][sfb][i] )
            {
                over++;
                l3_xmin->s[gr][ch][sfb][i] *= ifqstep2;
                scalefac->s[gr][ch][sfb][i]++;
#ifdef DEBUGSC
                printf( "cod_info->scalefac[%d][%d] = %d (amp_scale)\n",
                        i,sfb,scalefac->s[gr][ch][sfb][i] );
#endif
                start = scalefac_band_short[sfb];
                end   = scalefac_band_short[sfb+1];
                for ( l = start; l < end; l++ )
                    (*xr_s)[l][i] *= ifqstep;
            }
    return over;
}




/*************************************************************************/
/*            quantize                                                   */
/*************************************************************************/

/*
  Function: Quantization of the vector xr ( -> ix)
*/

void quantize( double xr[576], int ix[576], gr_info *cod_info )
{
    int i, b, l_end, s_start;
    double ostep, step, quantizerStepSize;
    static int inited = 0;

    D192_3 *xr_s;
    I192_3 *ix_s;

    if (!inited) {
      inited = 1;
      /* appr_pow_0_75_init(); */
      init_pow_nint();
    }
    
    xr_s = (D192_3 *) xr;
    ix_s = (I192_3 *) ix;

    quantizerStepSize = (double) cod_info->quantizerStepSize;

    for ( i = 0; i < 576; i++ )
	ix[i] = 0;

    if ( cod_info->quantizerStepSize == 0.0 )
	step = 1.0;
    else
	step = pow ( 2.0, quantizerStepSize * 0.25 );

    if ( cod_info->window_switching_flag != 0 && cod_info->block_type == 2 )
	if ( cod_info->mixed_block_flag == 0 )
	{
	    l_end = 0;
	    s_start = 0;
	}
	else
	{
	    l_end = 18 * 2;
	    s_start = 6 * 2;
	}
    else
    {
	l_end = 576;
	s_start = 192;
    }

    ostep = 1.0/step;
    for ( i = 0; i < l_end; i++ ) {
      /* ix[i] = nint( appr_pow_0_75 (fabs(xr[i]) * ostep) - 0.0946 ); */
      /*      ix[i] = nint( pow (fabs(xr[i]) * ostep, 0.75) - 0.0946 ); */
      
      ix[i] = pow_nint (fabs(xr[i]) * ostep);
      
	/*	if ((ix[i] > biggest) || (ix[i] < 0)) {
	    printf ("%d   \n", ix[i]);
	    biggest = ix[i];
	 }
	 */
    }
    
    if ( s_start < 192 )
	for ( b = 0; b < 3; b++ )
	{
	    step = pow( 2.0, (quantizerStepSize + 8.0 * (double) cod_info->subblock_gain[b]) * 0.25 );
	    ostep = 1.0/step;
	    for ( i = s_start; i < 192; i++ )
	      /*(*ix_s)[i][b] = nint( pow(fabs((*xr_s)[i][b]) * ostep, 0.75) - 0.0946 );*/
		(*ix_s)[i][b] = pow_nint(fabs((*xr_s)[i][b]) * ostep);
	}
}




/*************************************************************************/
/*            ix_max                                                     */
/*************************************************************************/

/*
  Function: Calculate the maximum of ix from 0 to 575
*/

int ix_max( int ix[576], unsigned int begin, unsigned int end )
{
    int i, max = 0;

    for ( i = begin; i < end; i++ )
    {
        int x = abs( ix[i] );
        if ( x > max )
            max = x;
    }
    return max;
}


/*************************************************************************/
/*            xr_max                                                     */
/*************************************************************************/

/*
  Function: Calculate the maximum of xr[576]  from 0 to 575
*/

double xr_max( double xr[576], unsigned int begin, unsigned int end )
{
    int i;
    double max = 0.0, temp;

    for ( i = begin; i < end; i++ )
        if( (temp = fabs(xr[i])) > max )
	    max = temp;
    return max;
}



/*        Noiseless coding  -- Huffman coding   */


/*************************************************************************/
/*            calc_runlen                                                */
/*************************************************************************/

/*
Function: Calculation of rzero, count1, big_values
(Partitions ix into big values, quadruples and zeros).
*/

void calc_runlen( int ix[576], gr_info *cod_info )
{
    int i;
    int rzero = 0; 

    if ( cod_info->window_switching_flag && (cod_info->block_type == 2) )
    {  /* short blocks */
        cod_info->count1 = 0;
        cod_info->big_values = 288;
    }
    else
    {
        for ( i = 576; i > 1; i -= 2 )
            if ( ix[i-1] == 0 && ix[i-2] == 0 )
                rzero++;
            else
                break;
        
        cod_info->count1 = 0 ;
        for ( ; i > 3; i -= 4 )
            if ( abs(ix[i-1]) <= 1
              && abs(ix[i-2]) <= 1
              && abs(ix[i-3]) <= 1
              && abs(ix[i-4]) <= 1 )
                cod_info->count1++;
            else
                break;
        
        cod_info->big_values = i/2;
    }
    assert( (2 * rzero + 4 * cod_info->count1 + 2 * cod_info->big_values) == 576 );
}



/*************************************************************************/
/*            count1_bitcount                                            */
/*************************************************************************/

/*
  Determines the number of bits to encode the quadruples.
*/

int count1_bitcount( int ix[ 576 ], gr_info *cod_info )
{
    int abs_and_sign( int *x );

    int p, i, k, bitsum_count1;
    int v, w, x, y, signbits;
    int sum0 = 0, sum1 = 0;

    for ( i = cod_info->big_values * 2, k = 0; k < cod_info->count1; i += 4, k++ )
    { 
      v = abs(ix[i]);
      w = abs(ix[i+1]);
      x = abs(ix[i+2]);
      y = abs(ix[i+3]);


      /*
        v = ix[ i ];
        w = ix[ i + 1 ];
        x = ix[ i + 2 ];
        y = ix[ i + 3 ];
        

        abs_and_sign( &v );
        abs_and_sign( &w );
        abs_and_sign( &x );
        abs_and_sign( &y );
	*/
        p = v + (w << 1) + (x << 2) + (y << 3);
        
        signbits = 0;

        if ( v != 0 )
            signbits ++;
        if ( w != 0 )
            signbits ++;
        if ( x != 0 )
            signbits ++;
        if ( y != 0 )
            signbits ++;

        sum0 += signbits;
        sum1 += signbits;

        sum0 += ht[ 32 ].hlen[ p ];
        sum1 += ht[ 33 ].hlen[ p ];
    }

    if ( sum0 < sum1 )
    {
        bitsum_count1 = sum0;
        cod_info->count1table_select = 0;
    }
    else
    {
        bitsum_count1 = sum1;
        cod_info->count1table_select = 1;
    }
    return( bitsum_count1 );
}





struct
{
    unsigned region0_count;
    unsigned region1_count;
} subdv_table[ 23 ] =
{
{0, 0}, /* 0 bands */
{0, 0}, /* 1 bands */
{0, 0}, /* 2 bands */
{0, 0}, /* 3 bands */
{0, 0}, /* 4 bands */
{0, 1}, /* 5 bands */
{1, 1}, /* 6 bands */
{1, 1}, /* 7 bands */
{1, 2}, /* 8 bands */
{2, 2}, /* 9 bands */
{2, 3}, /* 10 bands */
{2, 3}, /* 11 bands */
{3, 4}, /* 12 bands */
{3, 4}, /* 13 bands */
{3, 4}, /* 14 bands */
{4, 5}, /* 15 bands */
{4, 5}, /* 16 bands */
{4, 6}, /* 17 bands */
{5, 6}, /* 18 bands */
{5, 6}, /* 19 bands */
{5, 7}, /* 20 bands */
{6, 7}, /* 21 bands */
{6, 7}, /* 22 bands */
};




/*************************************************************************/
/*            subdivide                                                  */
/*************************************************************************/

/* presumable subdivides the bigvalue region which will
   use separate Huffman tables.
*/

void subdivide( gr_info *cod_info )
{
    int scfb_anz = 0;
    int bigvalues_region;
    
    if ( cod_info->big_values == 0 )
    { /* no big_values region */
        cod_info->region0_count = 0;
        cod_info->region1_count = 0;
    }
    else
    {
        bigvalues_region = 2 * cod_info->big_values;

        if ( (cod_info->window_switching_flag == 0) )
        { /* long blocks */
            int thiscount, index;
            /* Calculate scfb_anz */
            while ( scalefac_band_long[scfb_anz] < bigvalues_region )
                scfb_anz++;
            assert( scfb_anz < 23 );

            cod_info->region0_count = subdv_table[scfb_anz].region0_count;
            thiscount = cod_info->region0_count;
            index = thiscount + 1;
            while ( thiscount && (scalefac_band_long[index] > bigvalues_region) )
            {
                thiscount -= 1;
                index -= 1;
            }
            cod_info->region0_count = thiscount;

            cod_info->region1_count = subdv_table[scfb_anz].region1_count;
            index = cod_info->region0_count + cod_info->region1_count + 2;
            thiscount = cod_info->region1_count;
            while ( thiscount && (scalefac_band_long[index] > bigvalues_region) )
            {
                thiscount -= 1;
                index -= 1;
            }
            cod_info->region1_count = thiscount;
            cod_info->address1 = scalefac_band_long[cod_info->region0_count+1];
            cod_info->address2 = scalefac_band_long[cod_info->region0_count
                                                    + cod_info->region1_count + 2 ];
            cod_info->address3 = bigvalues_region;
        }
        else
        {
            if ( (cod_info->block_type == 2) && (cod_info->mixed_block_flag == 0) )
            { 
                cod_info->region0_count =  8;
                cod_info->region1_count =  36;
                cod_info->address1 = 36;
                cod_info->address2 = bigvalues_region;
                cod_info->address3 = 0;  
            }
            else
            {
                cod_info->region0_count = 7;
                cod_info->region1_count = 13;
                cod_info->address1 = scalefac_band_long[cod_info->region0_count+1];
                cod_info->address2 = bigvalues_region;
                cod_info->address3 = 0;
            }
        }
    }
}




/*************************************************************************/
/*            bigv_tab_select                                            */
/*************************************************************************/

/*
/*  Function: Select huffman code tables for bigvalues regions 
*/

void bigv_tab_select( int ix[576], gr_info *cod_info )
{
    /* int max; */

    cod_info->table_select[0] = 0;
    cod_info->table_select[1] = 0;
    cod_info->table_select[2] = 0;
    
    if ( cod_info->window_switching_flag && (cod_info->block_type == 2) )
    {
        /*
          Within each scalefactor band, data is given for successive
          time windows, beginning with window 0 and ending with window 2.
          Within each window, the quantized values are then arranged in
          order of increasing frequency...
          */
        int sfb, window, line, start, end, max1, max2, x, y;
        int region1Start;
        int *pmax;

        region1Start = 12;
        max1 = max2 = 0;

        for ( sfb = 0; sfb < 13; sfb++ )
        {
            start = scalefac_band_short[ sfb ];
            end   = scalefac_band_short[ sfb+1 ];
            
            if ( start < region1Start )
                pmax = &max1;
            else
                pmax = &max2;
            
            for ( window = 0; window < 3; window++ )
                for ( line = start; line < end; line += 2 )
                {
                    assert( line >= 0 );
                    assert( line < 576 );
                    x = abs( ix[ (line * 3) + window ] );
                    y = abs( ix[ ((line + 1) * 3) + window ]);
                    *pmax = *pmax > x ? *pmax : x;
                    *pmax = *pmax > y ? *pmax : y;
                }
        }
        cod_info->table_select[0] = choose_table( max1 );
        cod_info->table_select[1] = choose_table( max2 );
    }
    else
    {
        if ( cod_info->address1 > 0 )
            cod_info->table_select[0] = new_choose_table( ix, 0, cod_info->address1 );

        if ( cod_info->address2 > cod_info->address1 )
            cod_info->table_select[1] = new_choose_table( ix, cod_info->address1, cod_info->address2 );

        if ( cod_info->big_values * 2 > cod_info->address2 )
            cod_info->table_select[2] = new_choose_table( ix, cod_info->address2, cod_info->big_values * 2 );
    }
}




/*************************************************************************/
/*            new_choose table                                           */
/*************************************************************************/

/*
  Choose the Huffman table that will encode ix[begin..end] with
  the fewest bits.

  Note: This code contains knowledge about the sizes and characteristics
  of the Huffman tables as defined in the IS (Table B.7), and will not work
  with any arbitrary tables.
*/

int new_choose_table( int ix[576], unsigned int begin, unsigned int end )
{
    int i, max;
    int choice[ 2 ];
    int sum[ 2 ];

    max = ix_max( ix, begin, end );

    if ( max == 0 )
        return 0;
    
    max = abs( max );

    choice[ 0 ] = 0;
    choice[ 1 ] = 0;

    if ( max < 15 )
    {
	/* try tables with no linbits */
        for ( i = 0; i < 14; i++ )
            if ( ht[i].xlen > max )
	    {
		choice[ 0 ] = i;
                break;
	    }
	assert( choice[0] );

	sum[ 0 ] = count_bit( ix, begin, end, choice[0] );

	switch ( choice[0] )
	{
	  case 2:
	    sum[ 1 ] = count_bit( ix, begin, end, 3 );
	    if ( sum[1] <= sum[0] )
		choice[ 0 ] = 3;
	    break;

	  case 5:
	    sum[ 1 ] = count_bit( ix, begin, end, 6 );
	    if ( sum[1] <= sum[0] )
		choice[ 0 ] = 6;
	    break;

	  case 7:
	    sum[ 1 ] = count_bit( ix, begin, end, 8 );
	    if ( sum[1] <= sum[0] )
	    {
		choice[ 0 ] = 8;
		sum[ 0 ] = sum[ 1 ];
	    }
	    sum[ 1 ] = count_bit( ix, begin, end, 9 );
	    if ( sum[1] <= sum[0] )
		choice[ 0 ] = 9;
	    break;

	  case 10:
	    sum[ 1 ] = count_bit( ix, begin, end, 11 );
	    if ( sum[1] <= sum[0] )
	    {
		choice[ 0 ] = 11;
		sum[ 0 ] = sum[ 1 ];
	    }
	    sum[ 1 ] = count_bit( ix, begin, end, 12 );
	    if ( sum[1] <= sum[0] )
		choice[ 0 ] = 12;
	    break;

	  case 13:
	    sum[ 1 ] = count_bit( ix, begin, end, 15 );
	    if ( sum[1] <= sum[0] )
		choice[ 0 ] = 15;
	    break;

	  default:
	    break;
	}
    }
    else
    {
	/* try tables with linbits */
	max -= 15;
	
	for ( i = 15; i < 24; i++ )
	{
	    if ( ht[i].linmax >= max )
	    {
		choice[ 0 ] = i;
		break;
	    }
	}
	for ( i = 24; i < 32; i++ )
	{
	    if ( ht[i].linmax >= max )
	    {
		choice[ 1 ] = i;
		break;
	    }
	}
	assert( choice[0] );
	assert( choice[1] );
	
	sum[ 0 ] = count_bit( ix, begin, end, choice[0] );
	sum[ 1 ] = count_bit( ix, begin, end, choice[1] );
	if ( sum[1] < sum[0] )
	    choice[ 0 ] = choice[ 1 ];
    }
    return choice[ 0 ];
}



/*************************************************************************/
/*            choose table                                               */
/*************************************************************************/

int choose_table( int max )
{
    int  i, choice;

    if ( max == 0 )
        return 0;
    
    max = abs( max );    
    choice = 0;

    if ( max < 15 )
    {
        for ( i = 0; i < 15; i++ )
        {
            if ( ht[i].xlen > max )
            {
		choice = i;
		break;
            }
        }
    }
    else
    {	
	max -= 15;
	for (i = 15; i < 32; i++ )
	{
	    if ( ht[i].linmax >= max )
	    {
		choice = i;
		break;
	    }
	}
    }
    assert( choice );
    return choice;
}


/*************************************************************************/
/*            bigv_bitcount                                              */
/*************************************************************************/

/*
Function: Count the number of bits necessary to code the bigvalues region.
*/

int bigv_bitcount( int ix[576], gr_info *gi )
{
    int bits = 0;
    
    if ( gi->window_switching_flag && gi->block_type == 2 )
    {
        /*
          Within each scalefactor band, data is given for successive
          time windows, beginning with window 0 and ending with window 2.
          Within each window, the quantized values are then arranged in
          order of increasing frequency...
          */
        int sfb, window, line, start, end;
        I192_3 *ix_s;

        if ( gi->mixed_block_flag )
        {
            unsigned int table;

            if ( (table = gi->table_select[0]) != 0 )
                bits += count_bit( ix, 0, gi->address1, table );
            sfb = 2;
        }
        else
            sfb = 0;

        ix_s = (I192_3 *) &ix[0];

        for ( ; sfb < 13; sfb++ )
        {
            unsigned tableindex = 100;

            start = scalefac_band_short[ sfb ];
            end   = scalefac_band_short[ sfb+1 ];

            if ( start < 12 )
                tableindex = gi->table_select[ 0 ];
            else
                tableindex = gi->table_select[ 1 ];
            assert( tableindex < 32 );

            for ( window = 0; window < 3; window++ )
                for ( line = start; line < end; line += 2 )
                {
                    int x = (*ix_s)[line][window];
                    int y = (*ix_s)[line + 1][window];
                    bits += HuffmanCode( tableindex, x, y, NULL );
                }
        }
    }
    else
    {
        unsigned int table;
        
        if( (table=gi->table_select[0]) != 0 )  /* region0 */ 
            bits += count_bit(ix, 0, gi->address1, table );
        if( (table=gi->table_select[1]) != 0 )  /* region1 */ 
            bits += count_bit(ix, gi->address1, gi->address2, table );
        if( (table=gi->table_select[2]) != 0 )  /* region2 */ 
            bits += count_bit(ix, gi->address2, gi->address3, table );
    }
    return bits;
}



#ifndef HAVE_NINT
int
nint( double in )
{
    int    temp;

    if( in < 0 )  temp = (int)(in - 0.5);
    else    temp = (int)(in + 0.5);

    return(temp);
}

double
aint(double in) {
	return((long) in);
}
#endif



/*
  Seymour's comment:  Jan 8 1995
  When mixed_block_flag is set, the low subbands 0-1 undergo the long
  window transform and are each split into 18 frequency lines, while
  the remaining 30 subbands undergo the short window transform and are
  each split into 6 frequency lines. A problem now arises, as neither
  the short or long scale factor bands apply to this mixed spectrum.
  The standard resolves this situation by using the first 8 long scale
  factor bands for the low spectrum and the short scale factor bands
  in the range of 3 to 11 (inclusive) for the remaining frequency lines.
  These scale factor bands do not match exactly to the 0-1 subbands
  for all sampling frequencies (32,44.1 and 48 kHz); however they
  were designed so that there would not be a frequency gap or overlap
  at the switch over point. (Note multiply short frequency lines by 3
  to account for wider frequency line.) 
  */



/*************************************************************************/
/*            gr_deco                                                    */
/*************************************************************************/

void gr_deco( gr_info *cod_info )
{
    if ( cod_info->window_switching_flag != 0 && cod_info->block_type == 2 )
        if ( cod_info->mixed_block_flag == 0 )
        {
            cod_info->sfb_lmax = 0; /* No sb*/
            cod_info->sfb_smax = 0;
        }
        else
        {
            cod_info->sfb_lmax = 8;
            cod_info->sfb_smax = 3;
        }
    else
    {
        cod_info->sfb_lmax = SFB_LMAX - 1;
        cod_info->sfb_smax = SFB_SMAX - 1;    /* No sb */
    }
}








/* The following optional code written by Seymour Shlien
   will speed up the outer_loop code which is called
   by iteration_loop. When BIN_SEARCH is defined, the
   outer_loop function precedes the call to the function inner_loop
   with a call to bin_search gain defined below, which
   returns a good starting quantizerStepSize.
*/

#if defined(BIN_SEARCH) || defined(PERFORM) 
int count_bits(ix,cod_info)  
int  *ix; /*  I576  *ix; */
gr_info *cod_info;
{
int bits,max;
  calc_runlen(ix,cod_info);		/*rzero,count1,big_values*/
  max = ix_max( ix, 0,576);
  if(max > 8192) return 100000;         /* report unsuitable quantizer */
  bits = count1_bitcount(ix, cod_info); /*count1_table selection*/
  subdivide(cod_info);			/* bigvalues sfb division */
  bigv_tab_select(ix,cod_info);		/* codebook selection*/
  bits += bigv_bitcount(ix,cod_info);	/* bit count */
/* printf("\nglobal_gain = %f  bits= %d ",cod_info->quantizerStepSize,bits);*/
return bits;
}
#endif


#ifdef BIN_SEARCH

int bin_search_StepSize(int desired_rate, double start, int *ix,
           double xrs[576], gr_info * cod_info)
{
double top,bot,next,last;
int bit;
top = start;
bot = 200;
next = start;
do
  {
  last = next;
  next = aint((top+bot)/2.0);
  cod_info->quantizerStepSize = next;
  quantize(xrs,ix,cod_info);
  bit = count_bits(ix,cod_info);
  if (bit>desired_rate) top = next;
  else bot = next;
/*  printf("\n%f %f %f %d %d",next, top,bot,bit,desired_rate);*/
  }
  while ((bit != desired_rate) && fabs(last - next) > 1.0);
return next;
}

#endif




#ifdef PERFORM

/* The following code is used for exposing some problems with
   the outer_loop code. PERFORM should be defined to the
   frame number you wish to have additional output recorded
   in the file encode.log  - Seymour Shlien 14-Jan-97
*/

/*
float worst_xfsf_to_xmin_ratio(l3_xmin,xfsf,block_type,gr,ch)
double xfsf[4][CBLIMIT];
III_psy_xmin *l3_xmin;
int block_type,gr,ch;
*/
float worst_xfsf_to_xmin_ratio(III_psy_xmin *l3_xmin, double xfsf[4][CBLIMIT]
                             ,int block_type,int gr,int ch)

{
float ratio,maxratio;
int i,j;
maxratio =-100.0;
if (block_type != 2)
  for(i=0;i<21;i++)
    {
    ratio = 10.*log10(xfsf[0][i] /l3_xmin->l[gr][ch][i]);
    if (ratio > maxratio) maxratio = ratio;
    }
else
  {
  for(j=0;j<3;j++)
/* for(i = cod_info->sfb_smax; i <SFB_SMAX; i++) */
     for(i = 0; i <11; i++)
    {
    ratio = 10.*log10(xfsf[j+1][i] /l3_xmin->s[gr][ch][i][j]);
    if (ratio > maxratio) maxratio = ratio;
    }
  }
return maxratio;
}
 

print_ratios(handle_out,l3_xmin,xfsf,block_type,gr,ch)
FILE *handle_out;
double xfsf[4][CBLIMIT];
III_psy_xmin *l3_xmin;
int gr,ch;
int block_type;
{
float ratio;
int i,j;
if(block_type !=2)
 for (i=0;i<21;i++)
  {
   ratio = 100.0; /* signals undefined value in output */
   if(l3_xmin->l[gr][ch][i] >1.0e-20)
     ratio = 10.*log10(xfsf[0][i] /l3_xmin->l[gr][ch][i]);
   fprintf(handle_out,"%6.2f ",ratio);
   if(i%5==4) fprintf(handle_out,"\n");
  }
else
  
  for(j=0;j<3;j++)
   { fprintf(handle_out,"\n  block %d\n",j);
     for(i = 0; i <11; i++)
      {
      ratio = 10.*log10(xfsf[j+1][i] /l3_xmin->s[gr][ch][i][j]);
      fprintf(handle_out,"%6.2f ",ratio);
      if(i%5==4) fprintf(handle_out,"\n");
      }
   }
fprintf(handle_out,"\n");
}


print_scalefacs(handle_out,scalefac,block_type,gr,ch)
FILE *handle_out;
III_scalefac_t *scalefac;
int gr,ch;
int block_type;
{
int sfb,j;

if(block_type !=2)
 for ( sfb = 0; sfb < 21; sfb++ )
   {
   fprintf(handle_out,"%6d ", scalefac->l[gr][ch][sfb]);    
   if(sfb%5==4) fprintf(handle_out,"\n");
   }
else
 for (j=0;j<3;j++)
  {
  fprintf(handle_out,"\n block %d\n",j);
  for (sfb=0;sfb<11;sfb++)
   {
   fprintf(handle_out,"%6d ",scalefac->s[gr][ch][sfb][j]);
   if(sfb%5==4) fprintf(handle_out,"\n");
   }
  }
fprintf(handle_out,"\n");
}


print_quantized_values(FILE *handle, int ix[576], gr_info *cod_info)
{
int sfb,start,end,i,bw;
for (sfb=0;sfb<cod_info->sfb_lmax;sfb++)
  {
  start = scalefac_band_long[sfb];
  end = scalefac_band_long[sfb+1];
  bw = end - start;
  fprintf(handle,"scalefac band %d from %d to %d\n",sfb,start,end);
  for (i=0;i<bw;i++)
    {
    fprintf(handle,"%8d",ix[start+i]);
    if(i%5==4) fprintf(handle,"\n");
    }
  fprintf(handle,"\n");
  }
}


test_inner_loop(double xr[2][2][576], int l3_enc[2][2][576], int max_bits,
            gr_info *cod_info, int gr, int ch, double xfsf[4][CBLIMIT], 
            III_psy_xmin  *l3_xmin) 
{
    int bits, c1bits, bvbits;
    double *xrs;  /*  D576  *xr; */
    int *ix;  /*  I576  *ix; */
    xrs = &xr[gr][ch][0];
    ix = l3_enc[gr][ch];
 
    assert( max_bits >= 0 );
    cod_info->quantizerStepSize -= 1.0;;
    do
    {
        do
        {
            cod_info->quantizerStepSize += 1.0;
            quantize( xrs, ix, cod_info );
        }  
        while ( ix_max(ix, 0, 576) > 8191 + 14 ); /* within table range? */
 
    bits = count_bits(ix,cod_info);       

        if(frameNum == PERFORM)
        {
          fprintf(log_output,"StepSize=%f bits = %d huff_bits = %d\n",
            cod_info->quantizerStepSize,bits,max_bits);
          calc_noise( &xr[gr][ch][0], &l3_enc[gr][ch][0], cod_info, xfsf );
                                             /* distortion calculation */
          print_ratios(log_output,l3_xmin,xfsf,cod_info->block_type,gr,ch);
          fprintf(log_output,"\n\n");
          }

    }
    while ( bits > max_bits );
 
    return bits;
}
 
#endif PERFORM 
