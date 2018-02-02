/**********************************************************************
 * ISO MPEG Audio Subgroup Software Simulation Group (1996)
 * ISO 13818-3 MPEG-2 Audio Encoder - Lower Sampling Frequency Extension
 *
 * $Id: l3bitstream.h,v 1.1 1996/02/14 04:04:23 rowlands Exp $
 *
 * $Log: l3bitstream.h,v $
 * Revision 1.1  1996/02/14 04:04:23  rowlands
 * Initial revision
 *
 * Received from Mike Coleman
 **********************************************************************/

#ifndef L3_BITSTREAM_H
#define L3_BITSTREAM_H

#include "common.h"
#include "encoder.h"
#include "formatBitstream.h"

void III_format_bitstream( int              bitsPerFrame,
			   frame_params     *in_fr_ps,
			   int              l3_enc[2][2][576],
                           III_side_info_t  *l3_side,
			   III_scalefac_t   *scalefac,
			   Bit_stream_struc *in_bs,
			   double           (*xr)[2][576],
			   char             *ancillary,
			   int              anc_bits );

//void III_FlushBitstream();

int abs_and_sign( int *x ); /* returns signx and changes *x to abs(*x) */


#endif
