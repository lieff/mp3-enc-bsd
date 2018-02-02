/**********************************************************************
 * ISO MPEG Audio Subgroup Software Simulation Group (1996)
 * ISO 13818-3 MPEG-2 Audio Encoder - Lower Sampling Frequency Extension
 *
 * $Id: reservoir.h,v 1.1 1996/02/14 04:04:23 rowlands Exp $
 *
 * $Log: reservoir.h,v $
 * Revision 1.1  1996/02/14 04:04:23  rowlands
 * Initial revision
 *
 * Received from Mike Coleman
 **********************************************************************/
/*
  Revision History:

  Date        Programmer                Comment
  ==========  ========================= ===============================
  1995/09/06  mc@fivebats.com           created

*/

#ifndef RESERVOIR_H
#define RESERVOIR_H

void ResvFrameBegin( frame_params *fr_ps, III_side_info_t *l3_side, int mean_bits, int frameLength );
int  ResvMaxBits( frame_params *fr_ps, III_side_info_t *l3_side, double *pe, int mean_bits );
void ResvAdjust( frame_params *fr_ps, gr_info *gi, III_side_info_t *l3_side, int mean_bits );
void ResvFrameEnd( frame_params *fr_ps, III_side_info_t *l3_side, int mean_bits );

#endif
