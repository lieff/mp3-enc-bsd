/**********************************************************************
 * ISO MPEG Audio Subgroup Software Simulation Group (1996)
 * ISO 13818-3 MPEG-2 Audio Encoder - Lower Sampling Frequency Extension
 *
 * $Id: common.h,v 1.1 1996/02/14 04:04:23 rowlands Exp $
 *
 * $Log: common.h,v $
 * Revision 1.1  1996/02/14 04:04:23  rowlands
 * Initial revision
 *
 * Received from Mike Coleman
 **********************************************************************/
/**********************************************************************
 *   date   programmers         comment                               *
 * 2/25/91  Doulas Wong,        start of version 1.0 records          *
 *          Davis Pan                                                 *
 * 5/10/91  W. Joseph Carter    Reorganized & renamed all ".h" files  *
 *                              into "common.h" and "encoder.h".      *
 *                              Ported to Macintosh and Unix.         *
 *                              Added additional type definitions for *
 *                              AIFF, double/SANE and "bitstream.c".  *
 *                              Added function prototypes for more    *
 *                              rigorous type checking.               *
 * 27jun91  dpwe (Aware)        Added "alloc_*" defs & prototypes     *
 *                              Defined new struct 'frame_params'.    *
 *                              Changed info.stereo to info.mode_ext  *
 *                              #define constants for mode types      *
 *                              Prototype arguments if PROTO_ARGS     *
 * 5/28/91  Earle Jennings      added MS_DOS definition               *
 *                              MsDos function prototype declarations *
 * 7/10/91  Earle Jennings      added FLOAT definition as double      *
 *10/ 3/91  Don H. Lee          implemented CRC-16 error protection   *
 * 2/11/92  W. Joseph Carter    Ported new code to Macintosh.  Most   *
 *                              important fixes involved changing     *
 *                              16-bit ints to long or unsigned in    *
 *                              bit alloc routines for quant of 65535 *
 *                              and passing proper function args.     *
 *                              Removed "Other Joint Stereo" option   *
 *                              and made bitrate be total channel     *
 *                              bitrate, irrespective of the mode.    *
 *                              Fixed many small bugs & reorganized.  *
 *                              Modified some function prototypes.    *
 *                              Changed BUFFER_SIZE back to 4096.     *
 * 7/27/92  Michael Li          (re-)Ported to MS-DOS                 *
 * 7/27/92  Masahiro Iwadare    Ported to Convex                      *
 * 8/07/92  mc@tv.tek.com                                             *
 * 8/10/92  Amit Gulati         Ported to the AIX Platform (RS6000)   *
 *                              AIFF string constants redefined       *
 * 8/27/93 Seymour Shlien,      Fixes in Unix and MSDOS ports,        *
 *         Daniel Lauzon, and                                         *
 *         Bill Truerniet                                             *
 *--------------------------------------------------------------------*
 * 4/23/92  J. Pineda           Added code for Layer III.             *
 * 11/9/92  Amit Gulati         Added defines for layerIII stereo     *
 *                              modes.                                *
 *  8/24/93 Masahiro Iwadare    Included IS modification in Layer III.*
 *                              Changed for 1 pass decoding.          *
 *  9/07/93 Toshiyuki Ishino    Integrated Layer III with Ver 3.9.    *
 *--------------------------------------------------------------------*
 * 11/20/93 Masahiro Iwadare    Integrated Layer III with Ver 4.0.    *
 *--------------------------------------------------------------------*
 *  7/14/94 Juergen Koller      Fix for HPUX an IRIX in AIFF-Strings  *
 *  6/12/95 Soeren H. Nielsen   Bug fix in new_ext().                 *
 *  7/11/95 Soeren H. Nielsen   Added defs. for MPEG-2 LSF            *
 *--------------------------------------------------------------------*
 * 8/02/95  mc@fivebats.com     Changes to AIFF stuff for portability *
 *                              to little-endian machines             *
 **********************************************************************/
#ifndef COMMON_DOT_H
#define COMMON_DOT_H
/***********************************************************************
*
*  Global Conditional Compile Switches
*
***********************************************************************/

#if 0
#define      UNIX            /* Unix conditional compile switch */
#define      MACINTOSH       /* Macintosh conditional compile switch */
#define      MS_DOS          /* IBM PC conditional compile switch */
#define      MSC60           /* Compiled for MS_DOS with MSC v6.0 */
#define      AIX             /* AIX conditional compile switch    */
#define      CONVEX          /* CONVEX conditional compile switch */
#endif

#if defined(MSC60) 
#ifndef MS_DOS
#define MS_DOS
#endif
#ifndef PROTO_ARGS
#define PROTO_ARGS
#endif
#endif

#ifdef  UNIX
#define         TABLES_PATH     "tables"  /* to find data files */
/* name of environment variable holding path of table files */
#define         MPEGTABENV      "MPEGTABLES"
#define         PATH_SEPARATOR  "/"        /* how to build paths */
#endif  /* UNIX */

#ifdef  MACINTOSH
#define      TABLES_PATH ":tables:"  /* where to find data files */
#endif  /* MACINTOSH */

/* 
 * Don't define FAR to far unless you're willing to clean up the 
 * prototypes
 */
#define FAR /*far*/

#ifdef __STDC__
#ifndef PROTO_ARGS
#define PROTO_ARGS
#endif
#endif

#ifdef CONVEX
#define SEEK_SET        0
#define SEEK_CUR        1
#define SEEK_END        2
#endif

/* MS_DOS and VMS do not define TABLES_PATH, so OpenTableFile will default
   to finding the data files in the default directory */

/***********************************************************************
*
*  Global Include Files
*
***********************************************************************/

#include        <stdio.h>
#include        <string.h>
#include        <math.h>
#include "portableio.h"
#include "ieeefloat.h"

#ifdef  UNIX
#include        <unistd.h>
#endif  /* UNIX */

#ifdef __sgi
#include	<stdlib.h>
#endif

#ifdef  MACINTOSH
#include        <stdlib.h>
#include        <console.h>
#endif  /* MACINTOSH */

#ifdef  MS_DOS
#include        <stdlib.h>
#ifdef MSC60
#include        <memory.h>
#else
#include        <alloc.h>
#include        <mem.h>
#endif  /* MSC60 */
#endif  /* MS_DOS */

#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#endif

#ifndef EXIT_FAILURE
#define EXIT_FAILURE 1
#endif

/***********************************************************************
*
*  Global Definitions
*
***********************************************************************/

/* General Definitions */

#ifdef UNIX
#include <unistd.h>
#include <limits.h>
#endif /* UNIX */

#ifdef  MS_DOS
#define         FLOAT                   double
#else
#define         FLOAT                   float
#endif
#ifndef FALSE
#define         FALSE                   0
#endif

#ifndef TRUE
#define         TRUE                    1
#endif

#define         NULL_CHAR               '\0'

#define         MAX_U_32_NUM            0xFFFFFFFF
#ifndef PI
#define         PI                      3.14159265358979
#endif
#define         PI4                     PI/4
#define         PI64                    PI/64
#define         LN_TO_LOG10             0.2302585093

#define         VOL_REF_NUM             0
#define         MPEG_AUDIO_ID           1
#define		MPEG_PHASE2_LSF		0	/* 1995-07-11 SHN */
#define         MAC_WINDOW_SIZE         24

#define         MONO                    1
#define         STEREO                  2
#define         BITS_IN_A_BYTE          8
#define         WORD                    16

#ifdef PATH_MAX
#define         MAX_NAME_SIZE           (PATH_MAX+1)
#else
#define         MAX_NAME_SIZE           (256+1)
#endif

#define         SBLIMIT                 32
#define         SSLIMIT                 18
#define         FFT_SIZE                1024
#define         HAN_SIZE                512
#define         SCALE_BLOCK             12
#define         SCALE_RANGE             64
#define         SCALE                   32768
#define         CRC16_POLYNOMIAL        0x8005

/* MPEG Header Definitions - Mode Values */

#define         MPG_MD_STEREO           0
#define         MPG_MD_JOINT_STEREO     1
#define         MPG_MD_DUAL_CHANNEL     2
#define         MPG_MD_MONO             3

/* Mode Extention */

#define         MPG_MD_LR_LR             0
#define         MPG_MD_LR_I              1
#define         MPG_MD_MS_LR             2
#define         MPG_MD_MS_I              3


/* "bit_stream.h" Definitions */

#define         MINIMUM         4    /* Minimum size of the buffer in bytes */
#define         MAX_LENGTH      32   /* Maximum length of word written or
                                        read from bit stream */
#define         READ_MODE       0
#define         WRITE_MODE      1
#define         ALIGNING        8
#define         BINARY          0
#define         ASCII           1

#ifndef BS_FORMAT
#define         BS_FORMAT       ASCII /* BINARY or ASCII = 2x bytes */
#endif

#define         BUFFER_SIZE     4096

#define         MIN(A, B)       ((A) < (B) ? (A) : (B))
#define         MAX(A, B)       ((A) > (B) ? (A) : (B))

/***********************************************************************
*
*  Global Type Definitions
*
***********************************************************************/

/* Structure for Reading Layer II Allocation Tables from File */

typedef struct {
    unsigned int    steps;
    unsigned int    bits;
    unsigned int    group;
    unsigned int    quant;
} sb_alloc, *alloc_ptr;

typedef sb_alloc        al_table[SBLIMIT][16];

/* Header Information Structure */

typedef struct {
    int version;
    int lay;
    int error_protection;
    int bitrate_index;
    int sampling_frequency;
    int padding;
    int extension;
    int mode;
    int mode_ext;
    int copyright;
    int original;
    int emphasis;
} layer, *the_layer;

/* Parent Structure Interpreting some Frame Parameters in Header */

typedef struct {
    layer       *header;        /* raw header information */
    int         actual_mode;    /* when writing IS, may forget if 0 chs */
    al_table    *alloc;         /* bit allocation table read in */
    int         tab_num;        /* number of table as loaded */
    int         stereo;         /* 1 for mono, 2 for stereo */
    int         jsbound;        /* first band of joint stereo coding */
    int         sblimit;        /* total number of sub bands */
} frame_params;


/* AIFF Definitions */

#define IFF_ID_FORM 0x464f524d /* "FORM" */
#define IFF_ID_AIFF 0x41494646 /* "AIFF" */
#define IFF_ID_COMM 0x434f4d4d /* "COMM" */
#define IFF_ID_SSND 0x53534e44 /* "SSND" */
#define IFF_ID_MPEG 0x4d504547 /* "MPEG" */

#define AIFF_FORM_HEADER_SIZE 12
#define AIFF_SSND_HEADER_SIZE 16


typedef struct  blockAlign_struct {
    unsigned long   offset;
    unsigned long   blockSize;
} blockAlign;

typedef struct  IFF_AIFF_struct {
    short           numChannels;
    unsigned long   numSampleFrames;
    short           sampleSize;
    double          sampleRate;
    unsigned long   sampleType;
    blockAlign      blkAlgn;
} IFF_AIFF;

enum byte_order { order_unknown, order_bigEndian, order_littleEndian };
extern enum byte_order NativeByteOrder;

/* "bit_stream.h" Type Definitions */

typedef struct  bit_stream_struc {
    FILE        *pt;            /* pointer to bit stream device */
    unsigned char *buf;         /* bit stream buffer */
    int         buf_size;       /* size of buffer (in number of bytes) */
    long        totbit;         /* bit counter of bit stream */
    int         buf_byte_idx;   /* pointer to top byte in buffer */
    int         buf_bit_idx;    /* pointer to top bit of top byte in buffer */
    int         mode;           /* bit stream open in read or write mode */
    int         eob;            /* end of buffer index */
    int         eobs;           /* end of bit stream flag */
    char        format;
    
    /* format of file in rd mode (BINARY/ASCII) */
} Bit_stream_struc;

#include "l3side.h"

/***********************************************************************
*
*  Global Variable External Declarations
*
***********************************************************************/

extern char     *mode_names[4];
extern char     *layer_names[3];
extern char	*version_names[2];
extern double   s_freq[2][4];
extern int      bitrate[2][3][15];
extern double FAR multiple[64];

/***********************************************************************
*
*  Global Function Prototype Declarations
*
***********************************************************************/

/* The following functions are in the file "common.c" */

#ifdef  PROTO_ARGS
extern FILE           *OpenTableFile(char*);
extern int            read_bit_alloc(int, al_table*);
extern int            pick_table(frame_params*);
extern int            js_bound(int, int);
extern void           hdr_to_frps(frame_params*);
extern void           WriteHdr(frame_params*, FILE*);
extern void           WriteBitAlloc(unsigned int[2][SBLIMIT], frame_params*,
                        FILE*);
extern void           WriteScale(unsigned int[2][SBLIMIT],
                        unsigned int[2][SBLIMIT], unsigned int[2][3][SBLIMIT],
                        frame_params*, FILE*);
extern void           WriteSamples(int, unsigned int FAR [SBLIMIT],
                        unsigned int[SBLIMIT], frame_params*, FILE*);
extern int            NumericQ(char*);
extern int            BitrateIndex(int, int, int);
extern int            SmpFrqIndex(long, int*);
extern int            memcheck(char*, int, int);
extern void           FAR *mem_alloc(unsigned long, char*);
extern void           mem_free(void**);
extern void           double_to_extended(double*, char[10]);
extern void           extended_to_double(char[10], double*);
extern int            aiff_read_headers(FILE*, IFF_AIFF*);
extern int            aiff_seek_to_sound_data(FILE*);
extern int            aiff_write_headers(FILE*, IFF_AIFF*);
extern int            refill_buffer(Bit_stream_struc*);
extern void           empty_buffer(Bit_stream_struc*, int);
extern void           open_bit_stream_w(Bit_stream_struc*, char*, int);
extern void           open_bit_stream_r(Bit_stream_struc*, char*, int);
extern void           close_bit_stream_r(Bit_stream_struc*);
extern void           close_bit_stream_w(Bit_stream_struc*);
extern void           alloc_buffer(Bit_stream_struc*, int);
extern void           desalloc_buffer(Bit_stream_struc*);
extern void           back_track_buffer(Bit_stream_struc*, int);
extern unsigned int   get1bit(Bit_stream_struc*);
extern void           put1bit(Bit_stream_struc*, int);
extern unsigned long  look_ahead(Bit_stream_struc*, int);
extern unsigned long  getbits(Bit_stream_struc*, int);
extern void           putbits(Bit_stream_struc*, unsigned int, int);
extern void           byte_ali_putbits(Bit_stream_struc*, unsigned int, int);
extern unsigned long  byte_ali_getbits(Bit_stream_struc*, int);
extern unsigned long  sstell(Bit_stream_struc*);
extern int            end_bs(Bit_stream_struc*);
extern int            seek_sync(Bit_stream_struc*, long, int);
extern void           I_CRC_calc(frame_params*, unsigned int[2][SBLIMIT],
                        unsigned int*);
extern void           II_CRC_calc(frame_params*, unsigned int[2][SBLIMIT],
                        unsigned int[2][SBLIMIT], unsigned int*);
extern void           update_CRC(unsigned int, unsigned int, unsigned int*);
extern void           read_absthr(FLOAT*, int);
extern unsigned int   hget1bit(); /* MI */
extern unsigned long  hgetbits(int);
extern unsigned long  hsstell();
extern void           hputbuf(unsigned int,int);
extern enum byte_order DetermineByteOrder();
extern void SwapBytesInWords( short *loc, int words );


#ifdef  MACINTOSH
extern void           set_mac_file_attr(char[MAX_NAME_SIZE], short, OSType,
                        OSType);
#endif
#ifdef MS_DOS
extern void           new_ext(char *filename, char *extname, char *newname); 
#endif

#else
extern FILE           *OpenTableFile();
extern int            read_bit_alloc();
extern int            pick_table();
extern int            js_bound();
extern void           hdr_to_frps();
extern void           WriteHdr();
extern void           WriteBitAlloc();
extern void           WriteScale();
extern void           WriteSamples();
extern int            NumericQ();
extern int            BitrateIndex();
extern int            SmpFrqIndex();
extern int            memcheck();
extern void           FAR *mem_alloc();
extern void           mem_free();
extern void           double_to_extended();
extern void           extended_to_double();
extern int            aiff_read_headers();
extern int            aiff_seek_to_sound_data();
extern int            aiff_write_headers();
extern int            refill_buffer();
extern void           empty_buffer();
extern void           open_bit_stream_w();
extern void           open_bit_stream_r();
extern void           close_bit_stream_r();
extern void           close_bit_stream_w();
extern void           alloc_buffer();
extern void           desalloc_buffer();
extern void           back_track_buffer();
extern unsigned int   get1bit();
extern void           put1bit();
extern unsigned long  look_ahead();
extern unsigned long  getbits();
extern void           putbits();
extern void           byte_ali_putbits();
extern unsigned long  byte_ali_getbits();
extern unsigned long  sstell();
extern int            end_bs();
extern int            seek_sync();
extern void           I_CRC_calc();
extern void           II_CRC_calc();
extern void           update_CRC();
extern void           read_absthr();

extern unsigned int   hget1bit();
extern unsigned long  hgetbits();
extern unsigned long  hsstell();
extern void           hputbuf();

#ifdef MS_DOS
extern void           new_ext(); 
#endif
#endif
#endif
