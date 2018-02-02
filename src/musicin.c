/**********************************************************************
 * ISO MPEG Audio Subgroup Software Simulation Group (1996)
 * ISO 13818-3 MPEG-2 Audio Encoder - Lower Sampling Frequency Extension
 *
 * $Id: musicin.c,v 1.2 1997/01/19 22:28:29 rowlands Exp $
 *
 * $Log: musicin.c,v $
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
 * 3/01/91  Douglas Wong,       start of version 1.1 records          *
 *          Davis Pan                                                 *
 * 3/06/91  Douglas Wong,       rename: setup.h to endef.h            *
 *                              removed extraneous variables          *
 * 3/21/91  J.Georges Fritsch   introduction of the bit-stream        *
 *                              package. This package allows you      *
 *                              to generate the bit-stream in a       *
 *                              binary or ascii format                *
 * 3/31/91  Bill Aspromonte     replaced the read of the SB matrix    *
 *                              by an "code generated" one            *
 * 5/10/91  W. Joseph Carter    Ported to Macintosh and Unix.         *
 *                              Incorporated Jean-Georges Fritsch's   *
 *                              "bitstream.c" package.                *
 *                              Modified to strictly adhere to        *
 *                              encoded bitstream specs, including    *
 *                              "Berlin changes".                     *
 *                              Modified user interface dialog & code *
 *                              to accept any input & output          *
 *                              filenames desired.  Also added        *
 *                              de-emphasis prompt and final bail-out *
 *                              opportunity before encoding.          *
 *                              Added AIFF PCM sound file reading     *
 *                              capability.                           *
 *                              Modified PCM sound file handling to   *
 *                              process all incoming samples and fill *
 *                              out last encoded frame with zeros     *
 *                              (silence) if needed.                  *
 *                              Located and fixed numerous software   *
 *                              bugs and table data errors.           *
 * 27jun91  dpwe (Aware Inc)    Used new frame_params struct.         *
 *                              Clear all automatic arrays.           *
 *                              Changed some variable names,          *
 *                              simplified some code.                 *
 *                              Track number of bits actually sent.   *
 *                              Fixed padding slot, stereo bitrate    *
 *                              Added joint-stereo : scales L+R.      *
 * 6/12/91  Earle Jennings      added fix for MS_DOS in obtain_param  *
 * 6/13/91  Earle Jennings      added stack length adjustment before  *
 *                              main for MS_DOS                       *
 * 7/10/91  Earle Jennings      conversion of all float to FLOAT      *
 *                              port to MsDos from MacIntosh completed*
 * 8/ 8/91  Jens Spille         Change for MS-C6.00                   *
 * 8/22/91  Jens Spille         new obtain_parameters()               *
 *10/ 1/91  S.I. Sudharsanan,   Ported to IBM AIX platform.           *
 *          Don H. Lee,                                               *
 *          Peter W. Farrett                                          *
 *10/ 3/91  Don H. Lee          implemented CRC-16 error protection   *
 *                              newly introduced functions are        *
 *                              I_CRC_calc, II_CRC_calc and encode_CRC*
 *                              Additions and revisions are marked    *
 *                              with "dhl" for clarity                *
 *11/11/91 Katherine Wang       Documentation of code.                *
 *                                (variables in documentation are     *
 *                                surround by the # symbol, and an '*'*
 *                                denotes layer I or II versions)     *
 * 2/11/92  W. Joseph Carter    Ported new code to Macintosh.  Most   *
 *                              important fixes involved changing     *
 *                              16-bit ints to long or unsigned in    *
 *                              bit alloc routines for quant of 65535 *
 *                              and passing proper function args.     *
 *                              Removed "Other Joint Stereo" option   *
 *                              and made bitrate be total channel     *
 *                              bitrate, irrespective of the mode.    *
 *                              Fixed many small bugs & reorganized.  *
 * 2/25/92  Masahiro Iwadare    made code cleaner and more consistent *
 * 8/07/92  Mike Coleman        make exit() codes return error status *
 *                              made slight changes for portability   *
 *19 aug 92 Soren H. Nielsen    Changed MS-DOS file name extensions.  *
 * 8/25/92  Shaun Astarabadi    Replaced rint() function with explicit*
 *                              rounding for portability with MSDOS.  *
 * 9/22/92  jddevine@aware.com  Fixed _scale_factor_calc() calls.     *
 *10/19/92  Masahiro Iwadare    added info->mode and info->mode_ext   *
 *                              updates for AIFF format files         *
 * 3/10/93  Kevin Peterson      In parse_args, only set non default   *
 *                              bit rate if specified in arg list.    *
 *                              Use return value from aiff_read_hdrs  *
 *                              to fseek to start of sound data       *
 * 7/26/93  Davis Pan           fixed bug in printing info->mode_ext  *
 *                              value for joint stereo condition      *
 * 8/27/93 Seymour Shlien,      Fixes in Unix and MSDOS ports,        *
 *         Daniel Lauzon, and                                         *
 *         Bill Truerniet                                             *
 * 11/7/95 Soeren H. Nielsen    LSF added. Bug fix in MSDOS ext.      *
 * 8/02/95 mc@fivebats.com      Changed default bitrate selection so  *
 *                              it works with the new LSF stuff       *
 *10/01/95 mc@fivebats.com      Added layer3                          *
 **********************************************************************/

#ifdef MS_DOS
#include <dos.h>
#endif
#include <stdlib.h>
#include <sys/stat.h>
#include "common.h"
#include "encoder.h"
#include "l3psy.h"
#include "mdct.h"
#include "loop.h"
#include "l3bitstream.h"
#include <assert.h>

/* Global variable definitions for "musicin.c" */

FILE               *musicin;
Bit_stream_struc   bs;
char               *programName;
int		   iswav = 0, littleData = 0, streaming_input = 0;

/* Implementations */

/************************************************************************
*
* parse_args
*
* PURPOSE:  Sets encoding parameters to the specifications of the
* command line.  Default settings are used for parameters
* not specified in the command line.
*
* SEMANTICS:  The command line is parsed according to the following
* syntax:
*
* -l  is followed by the layer number
* -m  is followed by the mode
* -p  is followed by the psychoacoustic model number
* -s  is followed by the sampling rate
* -b  is followed by the total bitrate, irrespective of the mode
* -d  is followed by the emphasis flag
* -c  is followed by the copyright/no_copyright flag
* -o  is followed by the original/not_original flag
* -e  is followed by the error_protection on/off flag
* -L  is followed by the raw little endian flag
*
* If the input file is in AIFF format, the sampling frequency is read
* from the AIFF header.
*
* The input and output filenames are read into #inpath# and #outpath#.
*
************************************************************************/
 
void
parse_args(argc, argv, fr_ps, psy, num_samples, inPath, outPath)
int     argc;
char    **argv;
frame_params  *fr_ps;
int     *psy;
unsigned long *num_samples;
char    inPath[MAX_NAME_SIZE];
char    outPath[MAX_NAME_SIZE];
{
   FLOAT srate;
   int   brate;
   layer *info = fr_ps->header;
   int   err = 0, i = 0;
   IFF_AIFF pcm_aiff_data;
   long samplerate;
   long soundPosition;
   char ckwav[5];
 
   /* preset defaults */
   inPath[0] = '\0';   outPath[0] = '\0';
   info->lay = DFLT_LAY;
   switch(DFLT_MOD) {
      case 's': info->mode = MPG_MD_STEREO; info->mode_ext = 0; break;
      case 'd': info->mode = MPG_MD_DUAL_CHANNEL; info->mode_ext=0; break;
      case 'j': info->mode = MPG_MD_JOINT_STEREO; break;
      case 'm': info->mode = MPG_MD_MONO; info->mode_ext = 0; break;
      default:
         fprintf(stderr, "%s: Bad mode dflt %c\n", programName, DFLT_MOD);
         abort();
   }
   *psy = DFLT_PSY;
   if((info->sampling_frequency = SmpFrqIndex((long)(1000*DFLT_SFQ), &info->version)) < 0) {
      fprintf(stderr, "%s: bad sfrq default %.2f\n", programName, DFLT_SFQ);
      abort();
   }
  info->bitrate_index = 9;
  brate = 0;
   switch(DFLT_EMP) {
      case 'n': info->emphasis = 0; break;
      case '5': info->emphasis = 1; break;
      case 'c': info->emphasis = 3; break;
      default: 
         fprintf(stderr, "%s: Bad emph dflt %c\n", programName, DFLT_EMP);
         abort();
   }
   info->copyright = 0; info->original = 0; info->error_protection = FALSE;
 
   /* process args */
   while(++i<argc && err == 0) {
      char c, *token, *arg, *nextArg;
      int  argUsed;
 
      token = argv[i];
      if(*token++ == '-' && *token != '\0') {
         if(i+1 < argc) nextArg = argv[i+1];
         else           nextArg = "";
         argUsed = 0;
         while( (c = *token++) ) {
            if(*token /* NumericQ(token) */) arg = token;
            else                             arg = nextArg;
            switch(c) {
               case 'l':        info->lay = atoi(arg); argUsed = 1;
                  if(info->lay<1 || info->lay>3) {
                     fprintf(stderr,"%s: -l layer must be 1, 2, or 3, not %s\n",
                          programName, arg);
                     err = 1;
                  }
                  break;
               case 'm':        argUsed = 1;
                  if (*arg == 's')
                    { info->mode = MPG_MD_STEREO; info->mode_ext = 0; }
                  else if (*arg == 'd')
                    { info->mode = MPG_MD_DUAL_CHANNEL; info->mode_ext=0; }
                  else if (*arg == 'j')
                    { info->mode = MPG_MD_JOINT_STEREO; }
                  else if (*arg == 'm')
                    { info->mode = MPG_MD_MONO; info->mode_ext = 0; }
                  else {
                    fprintf(stderr,"%s: -m mode must be s/d/j/m not %s\n",
                            programName, arg);
                    err = 1;
                  }
                  break;
               case 'p':        *psy = atoi(arg); argUsed = 1;
                  if(*psy<1 || *psy>2) {
                     fprintf(stderr,"%s: -p model must be 1 or 2, not %s\n",
                             programName, arg);
                     err = 1;
                  }
                  break;

               case 's':
                  argUsed = 1;
                  srate = atof( arg );
                  /* samplerate = rint( 1000.0 * srate ); $A  */
                  samplerate = (long) (( 1000.0 * srate ) + 0.5);
                  if( (info->sampling_frequency =
          SmpFrqIndex((long) samplerate, &info->version)) < 0 )
                      err = 1;
                  break;
                  
               case 'b':        
        argUsed = 1;
      brate = atoi(arg); 
      break;
               case 'd':        argUsed = 1;
                  if (*arg == 'n')                    info->emphasis = 0;
                  else if (*arg == '5')               info->emphasis = 1;
                  else if (*arg == 'c')               info->emphasis = 3;
                  else {
                     fprintf(stderr,"%s: -d emp must be n/5/c not %s\n",
                             programName, arg);
                     err = 1;
                  }
                  break;
                case 'c':       info->copyright = 1; break;
                case 'o':       info->original  = 1; break;
                case 'e':       info->error_protection = TRUE; break;
	        case 'L':       littleData = TRUE; break;
                default:        fprintf(stderr,"%s: unrec option %c\n",
                                        programName, c);
                                err = 1; break;
            }
            if(argUsed) {
               if(arg == token)    token = "";   /* no more from token */
               else                ++i;          /* skip arg we used */
               arg = ""; argUsed = 0;
            }
         }
      }
      else {
         if(inPath[0] == '\0')       strcpy(inPath, argv[i]);
         else if(outPath[0] == '\0') strcpy(outPath, argv[i]);
         else {
            fprintf(stderr,"%s: excess arg %s\n", programName, argv[i]);
            err = 1;
         }
      }
   }

   if(err || inPath[0] == '\0') usage();  /* never returns */
 
   if(outPath[0] == '\0') {
#ifdef MS_DOS
      /* replace old extension with new one, 1992-08-19, 1995-06-12 shn */
      new_ext(inPath, DFLT_EXT, outPath);
#else
      strcpy(outPath, inPath);
      strcat(outPath, DFLT_EXT);
#endif
   }

   if (!strcmp (inPath, "-")) {
     musicin = stdin;
     streaming_input++;
   } else if ((musicin = fopen(inPath, "rb")) == NULL) {
      printf("Could not find \"%s\".\n", inPath);
      exit(1);
   }
 
   open_bit_stream_w(&bs, outPath, BUFFER_SIZE);

   if (!streaming_input) {
     if ((soundPosition = aiff_read_headers(musicin, &pcm_aiff_data)) != -1) {

        printf(">>> Using Audio IFF sound file headers\n");

        aiff_check(inPath, &pcm_aiff_data, &info->version);

        if (fseek(musicin, soundPosition, SEEK_SET) != 0) {
           printf("Could not seek to PCM sound data in \"%s\".\n", inPath);
           exit(1);
        }

        info->sampling_frequency = SmpFrqIndex((long)pcm_aiff_data.sampleRate, &info->version);
        printf(">>> %f Hz sampling frequency selected\n",
               pcm_aiff_data.sampleRate);

        /* Determine number of samples in sound file */
#ifndef MS_DOS
        *num_samples = pcm_aiff_data.numChannels *
                       pcm_aiff_data.numSampleFrames;
#else
        *num_samples = (long)(pcm_aiff_data.numChannels) *
                       (long)(pcm_aiff_data.numSampleFrames);
#endif
        if ( pcm_aiff_data.numChannels == 1 ) {
          info->mode = MPG_MD_MONO;
          info->mode_ext = 0;
        }
     }
     else {    /* Not using Audio IFF sound file headers. Assume PCM */
	struct stat sb;
	
        if (fseek(musicin, 8, SEEK_SET) != 0) {
           printf("Could not seek to PCM sound data in \"%s\".\n", inPath);
           exit(1);
        }

        fread(ckwav, 4, 1, musicin);
        if (!memcmp(ckwav, "WAVE", 4)) {
       	  iswav = 1;
          fseek(musicin, 0x2c, SEEK_SET); // skip .wav file header
        } else
          fseek(musicin, 0, SEEK_SET);

	if (fstat (fileno(musicin), &sb))
	  /* Declare sound file to have "infinite" number of samples. */
          *num_samples = MAX_U_32_NUM;
	else
	  *num_samples = (sb.st_size - (iswav ? 0x2c : 0)) / 2;
     }
   } else
     *num_samples = MAX_U_32_NUM;
   
   if ( brate == 0 )
    brate = bitrate[info->version][info->lay-1][9];
   if( (info->bitrate_index = BitrateIndex(info->lay, brate, info->version)) < 0) err=1;
   if(err || inPath[0] == '\0') usage();  /* never returns */

}

/************************************************************************
*
* print_config
*
* PURPOSE:  Prints the encoding parameters used
*
************************************************************************/
 
void print_config( frame_params *fr_ps, int *psy, char *inPath, char *outPath)
{
 layer *info = fr_ps->header;
 
   printf ("-------------------------------------------------------------------------------\n");
   printf("\tAlgorithm %s, layer %s/%s, bitrate %d kbps\n",
	  version_names[info->version],
	  layer_names[info->lay-1],
	  mode_names[info->mode],
	  bitrate[info->version][info->lay-1][info->bitrate_index]
	  );
   printf ("\tsampling frequency %.1f kHz,        psy. model %d\n",
	  s_freq[info->version][info->sampling_frequency], *psy);
   printf("\tinput  file: '%s'\n\toutput file: '%s'\n",
	  !strcmp (inPath, "-") ? "stdin" : inPath, outPath);
   printf ("-------------------------------------------------------------------------------\n");
}




 
/************************************************************************
*
* main
*
* PURPOSE:  MPEG I Encoder supporting layers 1 and 2, and 3, with
* psychoacoustic models 1 (MUSICAM) and 2 (AT&T)
*
* SEMANTICS:  One overlapping frame of audio of up to 2 channels are
* processed at a time in the following order:
* (associated routines are in parentheses)
*
* 1.  Filter sliding window of data to get 32 subband
* samples per channel.
* (window_subband,filter_subband)
*
* 2.  If joint stereo mode, combine left and right channels
* for subbands above #jsbound#.
* (*_combine_LR)
*
* 3.  Calculate scalefactors for the frame, and if layer 2,
* also calculate scalefactor select information.
* (*_scale_factor_calc)
*
* 4.  Calculate psychoacoustic masking levels using selected
* psychoacoustic model.
* (*_Psycho_One, psycho_anal)
*
* 5.  Perform iterative bit allocation for subbands with low
* mask_to_noise ratios using masking levels from step 4.
* (*_main_bit_allocation)
*
* 6.  If error protection flag is active, add redundancy for
* error protection.
* (*_CRC_calc)
*
* 7.  Pack bit allocation, scalefactors, and scalefactor select
* information (layer 2) onto bitstream.
* (*_encode_bit_alloc,*_encode_scale,II_transmission_pattern)
*
* 8.  Quantize subbands and pack them into bitstream
* (*_subband_quantization, *_sample_encoding)
*
************************************************************************/

int frameNum=0;

void main(argc, argv)
int     argc;
char    **argv;
{
    typedef double SBS[2][3][SCALE_BLOCK][SBLIMIT];
    SBS  FAR        *sb_sample;
    L3SBS  FAR        *l3_sb_sample;
    typedef double JSBS[3][SCALE_BLOCK][SBLIMIT];
    JSBS FAR        *j_sample;
    typedef double IN[2][HAN_SIZE];
    IN   FAR        *win_que;
    typedef unsigned int SUB[2][3][SCALE_BLOCK][SBLIMIT];
    SUB  FAR        *subband;
    
    frame_params fr_ps;
    layer info;
    char original_file_name[MAX_NAME_SIZE];
    char encoded_file_name[MAX_NAME_SIZE];
    short FAR **win_buf;
    static short FAR buffer[2][1152];
    static unsigned int bit_alloc[2][SBLIMIT], scfsi[2][SBLIMIT];
    static unsigned int scalar[2][3][SBLIMIT], j_scale[3][SBLIMIT];
    static double FAR ltmin[2][SBLIMIT], lgmin[2][SBLIMIT], max_sc[2][SBLIMIT];
    FLOAT snr32[32];
    short sam[2][1344]; /* was [1056]; */
    int whole_SpF, extra_slot = 0;
    double avg_slots_per_frame, frac_SpF, slot_lag;
    int model, stereo, error_protection;
    static unsigned int crc;
    int i, j, k, adb;
    unsigned long bitsPerSlot, samplesPerFrame;
    unsigned long frameBits, sentBits = 0;
    unsigned long num_samples, num_frames;
    
#ifdef  MACINTOSH
    argc = ccommand( &argv );
#endif
    

    /* Most large variables are declared dynamically to ensure
       compatibility with smaller machines */
    
    sb_sample = (SBS FAR *) mem_alloc(sizeof(SBS), "sb_sample");
    l3_sb_sample = (L3SBS FAR *) mem_alloc(sizeof(L3SBS), "l3_sb_sample");
    j_sample = (JSBS FAR *) mem_alloc(sizeof(JSBS), "j_sample");
    win_que = (IN FAR *) mem_alloc(sizeof(IN), "Win_que");
    subband = (SUB FAR *) mem_alloc(sizeof(SUB),"subband");
    win_buf = (short FAR **) mem_alloc(sizeof(short *)*2, "win_buf");
 
    /* clear buffers */
    memset((char *) buffer, 0, sizeof(buffer));
    memset((char *) bit_alloc, 0, sizeof(bit_alloc));
    memset((char *) scalar, 0, sizeof(scalar));
    memset((char *) j_scale, 0, sizeof(j_scale));
    memset((char *) scfsi, 0, sizeof(scfsi));
    memset((char *) ltmin, 0, sizeof(ltmin));
    memset((char *) lgmin, 0, sizeof(lgmin));
    memset((char *) max_sc, 0, sizeof(max_sc));
    memset((char *) snr32, 0, sizeof(snr32));
    memset((char *) sam, 0, sizeof(sam));
 
    fr_ps.header = &info;
    fr_ps.tab_num = -1;             /* no table loaded */
    fr_ps.alloc = NULL;
    info.version = MPEG_AUDIO_ID;   /* Default: MPEG-1 */

    programName = argv[0];
    if(argc==1) {     /* no command-line args */
      usage();
      exit(0);
    }
    else
	parse_args(argc, argv, &fr_ps, &model, &num_samples,
		   original_file_name, encoded_file_name);
    print_config(&fr_ps, &model,
                 original_file_name, encoded_file_name);
    
    hdr_to_frps(&fr_ps);
    stereo = fr_ps.stereo;
    error_protection = info.error_protection;
    
    if (info.lay == 1)
    { bitsPerSlot = 32; samplesPerFrame = 384;  }
    else 
	if ( info.lay == 2 )
	{ bitsPerSlot = 8;  samplesPerFrame = 1152; }
	else	
	{  /* layer 3 */
	    bitsPerSlot = 8;
	    samplesPerFrame = info.version == 1 ? 1152 : 576;
	    
	    /* Apologize for missing features */
	    if ( info.mode == MPG_MD_JOINT_STEREO )
	    {
		fprintf( stderr, "Sorry, joint stereo not yet available for layer3\n" );
		exit( 1 );
	    }
	    
	    if ( model != 2 )
	    {
		fprintf( stderr, "Sorry, psycho model 1 not available for layer3\n" );
		exit( 1 );
	    }
	}
    /* Figure average number of 'slots' per frame. */
    /* Bitrate means TOTAL for both channels, not per side. */
    avg_slots_per_frame = ((double)samplesPerFrame /
                           s_freq[info.version][info.sampling_frequency]) *
			   ((double)bitrate[info.version][info.lay-1][info.bitrate_index] /
			    (double)bitsPerSlot);
    whole_SpF = (int) avg_slots_per_frame;
    avg_slots_per_frame = whole_SpF;
#if 0
    printf("slots/frame = %d\n",whole_SpF);
#endif
    frac_SpF  = avg_slots_per_frame - (double)whole_SpF;
    slot_lag  = -frac_SpF;
#if 0
    printf("frac SpF=%.3f, tot bitrate=%d kbps, s freq=%.1f kHz\n",
           frac_SpF, bitrate[info.version][info.lay-1][info.bitrate_index],
           s_freq[info.version][info.sampling_frequency]);
#endif
    
    if (frac_SpF != 0)
	printf("Fractional number of slots, padding required\n");
    else info.padding = 0;

    num_frames = num_samples / samplesPerFrame / stereo;
    
    while ( get_audio(musicin, buffer, num_samples, stereo, &info) > 0 )
    {
        if (!(frameNum % 10))
	  if (num_samples == MAX_U_32_NUM) {
	    fprintf(stderr, "[%5lu]\r", frameNum);
	    fflush(stderr);
	  } else {
	    fprintf(stderr, "[%5lu/%5lu] (%02d%%)\r", frameNum, num_frames,
		    100*frameNum / num_frames);
	    fflush(stderr);
	  }
	frameNum++;
	
	win_buf[0] = &buffer[0][0];
	win_buf[1] = &buffer[1][0];
	if (frac_SpF != 0) {
	    if (slot_lag > (frac_SpF-1.0) ) {
		slot_lag -= frac_SpF;
		extra_slot = 0;
		info.padding = 0;
		/*  printf("No padding for this frame\n"); */
	    }
	    else {
		extra_slot = 1;
		info.padding = 1;
		slot_lag += (1-frac_SpF);
		/*  printf("Padding for this frame\n");    */
	    }
	}
	adb = (whole_SpF+extra_slot) * bitsPerSlot;
	
	switch (info.lay)
	{
	    
/***************************** Layer I **********************************/
	    
          case 1 :
	    for (j=0;j<SCALE_BLOCK;j++)
		for (k=0;k<stereo;k++) {
		    window_subband(&win_buf[k], &(*win_que)[k][0], k);
		    filter_subband(&(*win_que)[k][0], &(*sb_sample)[k][0][j][0]);
		}
	    
	    I_scale_factor_calc(*sb_sample, scalar, stereo);
	    if(fr_ps.actual_mode == MPG_MD_JOINT_STEREO) {
                I_combine_LR(*sb_sample, *j_sample);
                I_scale_factor_calc(j_sample, &j_scale, 1);
	    }
	    
	    put_scale(scalar, &fr_ps, max_sc);
	    
	    if (model == 1) I_Psycho_One(buffer, max_sc, ltmin, &fr_ps);
	    else {
                for (k=0;k<stereo;k++) {
		    psycho_anal(&buffer[k][0],&sam[k][0], k, info.lay, snr32,
				(FLOAT)s_freq[info.version][info.sampling_frequency]*1000);
		    for (i=0;i<SBLIMIT;i++) ltmin[k][i] = (double) snr32[i];
                }
	    }
	    
	    I_main_bit_allocation(ltmin, bit_alloc, &adb, &fr_ps);
	    
	    if (error_protection) I_CRC_calc(&fr_ps, bit_alloc, &crc);
	    
	    encode_info(&fr_ps, &bs);
	    
	    if (error_protection) encode_CRC(crc, &bs);
	    
	    I_encode_bit_alloc(bit_alloc, &fr_ps, &bs);
	    I_encode_scale(scalar, bit_alloc, &fr_ps, &bs);
	    I_subband_quantization(scalar, *sb_sample, j_scale, *j_sample,
				   bit_alloc, *subband, &fr_ps);
	    I_sample_encoding(*subband, bit_alloc, &fr_ps, &bs);
	    for (i=0;i<adb;i++) put1bit(&bs, 0);
	    break;
	    
/***************************** Layer 2 **********************************/
	    
          case 2 :
	    for (i=0;i<3;i++) for (j=0;j<SCALE_BLOCK;j++)
                for (k=0;k<stereo;k++) {
		    window_subband(&win_buf[k], &(*win_que)[k][0], k);
		    filter_subband(&(*win_que)[k][0], &(*sb_sample)[k][i][j][0]);
                }
	    
	    II_scale_factor_calc(*sb_sample, scalar, stereo, fr_ps.sblimit);
	    pick_scale(scalar, &fr_ps, max_sc);
	    if(fr_ps.actual_mode == MPG_MD_JOINT_STEREO) {
		II_combine_LR(*sb_sample, *j_sample, fr_ps.sblimit);
		II_scale_factor_calc(j_sample, &j_scale, 1, fr_ps.sblimit);
	    }       /* this way we calculate more mono than we need */
	    /* but it is cheap */
	    
	    if (model == 1) II_Psycho_One(buffer, max_sc, ltmin, &fr_ps);
	    else {
		for (k=0;k<stereo;k++) {
		    psycho_anal(&buffer[k][0],&sam[k][0], k, 
				info.lay, snr32,
				(FLOAT)s_freq[info.version][info.sampling_frequency]*1000);
		    for (i=0;i<SBLIMIT;i++) ltmin[k][i] = (double) snr32[i];
		}
	    }
	    
	    II_transmission_pattern(scalar, scfsi, &fr_ps);
	    II_main_bit_allocation(ltmin, scfsi, bit_alloc, &adb, &fr_ps);
	    
	    if (error_protection)
		II_CRC_calc(&fr_ps, bit_alloc, scfsi, &crc);
	    
	    encode_info(&fr_ps, &bs);
	    
	    if (error_protection) encode_CRC(crc, &bs);
	    
	    II_encode_bit_alloc(bit_alloc, &fr_ps, &bs);
	    II_encode_scale(bit_alloc, scfsi, scalar, &fr_ps, &bs);
	    II_subband_quantization(scalar, *sb_sample, j_scale,
				    *j_sample, bit_alloc, *subband, &fr_ps);
	    II_sample_encoding(*subband, bit_alloc, &fr_ps, &bs);
	    for (i=0;i<adb;i++) put1bit(&bs, 0);
	    break;
	    
/***************************** Layer 3 **********************************/

	  case 3:
	  {
	      /*
		large "auto" vars are static due to the Macintosh linker
	      */ 
	      static double xr[2][2][576];
	      static double xr_dec[2][2][576];
	      static double pe[2][2];
	      static int l3_enc[2][2][576];
	      static III_psy_ratio ratio;
	      static III_side_info_t l3_side;
	      static III_scalefac_t  scalefac;
	      int gr, mode_gr, ch;
	      int mean_bits, sideinfo_len;
	      
	      int bitsPerFrame = 8 * whole_SpF + (info.padding * 8);
	      mode_gr = (info.version == 1) ? 2 : 1;

	      /*
		determine the mean bitrate for main data
	      */
	      sideinfo_len = 32;
	      if ( info.version == 1 )
	      {   /* MPEG 1 */
		  if ( stereo == 1 )
		      sideinfo_len += 136;
		  else
		      sideinfo_len += 256;
	      }
	      else
	      {   /* MPEG 2 */
		  if ( stereo == 1 )
		      sideinfo_len += 72;
		  else
		      sideinfo_len += 136;
	      }
	      if ( info.error_protection )
		  sideinfo_len += 16;
	      mean_bits = (bitsPerFrame - sideinfo_len) / mode_gr;

	      /*
		psychoacoustic model
	      */
	      for ( gr = 0; gr < mode_gr; gr++ )
		  for ( ch = 0; ch < stereo; ch++ )
		  {
		      L3psycho_anal( &buffer[ch][gr*576], &sam[ch][0], ch, info.lay,
				     snr32, s_freq[info.version][info.sampling_frequency] * 1000.0,
				     &ratio.l[gr][ch][0], &ratio.s[gr][ch][0],
				     &pe[gr][ch], &l3_side.gr[gr].ch[ch].tt );
		  }

	      /*
		polyphase filtering
	      */
	      for( gr = 0; gr < mode_gr; gr++ )
		  for ( ch = 0; ch < stereo; ch++ )
		      for ( j = 0; j < 18; j++ )
		      {
			  window_subband( &win_buf[ch], &(*win_que)[ch][0], ch );
			  filter_subband( &(*win_que)[ch][0],  &(*l3_sb_sample)[ch][gr+1][j][0] );
		      }

	      /*
		apply mdct to the polyphase outputs
	      */
	      mdct_sub( l3_sb_sample, xr, stereo, &l3_side, mode_gr );
	      
	      /*
		bit and noise allocation
	      */
	      iteration_loop( pe, xr, &ratio, &l3_side, l3_enc, mean_bits,
			      stereo, xr_dec, &scalefac, &fr_ps, 0, bitsPerFrame );

	      /*
		write the frame to the bitstream
	      */
	      III_format_bitstream( bitsPerFrame, &fr_ps, l3_enc, &l3_side, &scalefac, &bs,
				    xr, NULL, 0 );
	  }
	    break;  /* end of layer 3 */
	    

	} /* end switch  */
	
	frameBits = sstell( &bs ) - sentBits;
	if ( frameBits % bitsPerSlot )   /* a program failure */
	    fprintf( stderr, "Sent %ld bits = %ld slots plus %ld\n",
		     frameBits, frameBits/bitsPerSlot,
		     frameBits%bitsPerSlot );
	sentBits += frameBits;

    }    

    if ( info.lay == 3 )
	III_FlushBitstream();

    close_bit_stream_w( &bs );

    printf("\nAvg slots/frame = %.3f; bits/sample = %.2f; bitrate = %.3f kbps\n",
           (FLOAT) sentBits / (frameNum * bitsPerSlot),
           (FLOAT) sentBits / (frameNum * samplesPerFrame),
           (FLOAT) sentBits / (frameNum * samplesPerFrame) *
           s_freq[info.version][info.sampling_frequency]);

    if (!streaming_input && fclose(musicin) != 0){
	printf("Could not close \"%s\".\n", original_file_name);
	exit(2);
    }

#ifdef  MACINTOSH
    set_mac_file_attr( encoded_file_name, VOL_REF_NUM, CREATOR_ENCODE,
		       FILETYPE_ENCODE );
#endif

    exit(0);
}
 
/************************************************************************
*
* usage
*
* PURPOSE:  Writes command line syntax to the file specified by #stderr#
*
************************************************************************/

void usage()  /* print syntax & exit */
{
    fprintf(stderr,
    "usage: %s [-l lay][-m mode][-p psy][-s sfrq][-b br][-d emp]\n",
            programName);
    fprintf(stderr,
    "          [-c][-o][-e] inputPCM [outBS]\n");
    fprintf(stderr,"where\n");
    fprintf(stderr," -l lay   use layer <lay> coding    (dflt %4u)\n",DFLT_LAY);
    fprintf(stderr," -m mode  channel mode : s/d/j/m    (dflt %4c)\n",DFLT_MOD);
    fprintf(stderr," -p psy   psychoacoustic model 1/2  (dflt %4u)\n",DFLT_PSY);
    fprintf(stderr," -s sfrq  input smpl rate in kHz    (dflt %4.1f)\n",DFLT_SFQ);
    fprintf(stderr," -b br    total bitrate in kbps     (dflt 128k)\n");
    fprintf(stderr," -d emp   de-emphasis n/5/c         (dflt %4c)\n",DFLT_EMP);
    fprintf(stderr," -L       PCM data is little endian (dflt big endian)\n");
    fprintf(stderr," -c       mark as copyright\n");
    fprintf(stderr," -o       mark as original\n");
    fprintf(stderr," -e       add error protection\n");
    fprintf(stderr," inputPCM input PCM sound file (WAV, standard or AIFF)\n");
    fprintf(stderr," outBS    output bit stream of encoded audio (dflt inName+%s)\n",
            DFLT_EXT);
    exit(1);
}

/************************************************************************
*
* aiff_check
*
* PURPOSE:  Checks AIFF header information to make sure it is valid.
*           Exits if not.
*
************************************************************************/

void aiff_check( char *file_name, IFF_AIFF *pcm_aiff_data, int *version)
{
    if (pcm_aiff_data->sampleType != IFF_ID_SSND) {
       printf("Sound data is not PCM in \"%s\".\n", file_name);
       exit(1);
    }

    if(SmpFrqIndex((long)pcm_aiff_data->sampleRate, version) < 0) {
       printf("in \"%s\".\n", file_name);
       exit(1);
    }

    if (pcm_aiff_data->sampleSize != sizeof(short) * BITS_IN_A_BYTE) {
        printf("Sound data is not %d bits in \"%s\".\n",
               sizeof(short) * BITS_IN_A_BYTE, file_name);
        exit(1);
    }

    if (pcm_aiff_data->numChannels != MONO &&
        pcm_aiff_data->numChannels != STEREO) {
       printf("Sound data is not mono or stereo in \"%s\".\n", file_name);
       exit(1);
    }

    if (pcm_aiff_data->blkAlgn.blockSize != 0) {
       printf("Block size is not %d bytes in \"%s\".\n", 0, file_name);
       exit(1);
    }

    if (pcm_aiff_data->blkAlgn.offset != 0) {
       printf("Block offset is not %d bytes in \"%s\".\n", 0, file_name);
       exit(1);
    }
}
