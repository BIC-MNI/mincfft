/* ----------------------------- MNI Header -----------------------------------
@NAME       : fft_support.h
@DESCRIPTION: defines and prototypes for fft support routines.
@CREATED    : Fri Nov  5 11:16:54 EST 1993 Louis Collins
@MODIFIED   : $Log: fft_support.h,v $
@MODIFIED   : Revision 1.1  2002-09-17 23:44:51  rotor
@MODIFIED   : ----------------------------------------------------------------------
@MODIFIED   : Initial entry of mincfft to CVS (to perhaps replace louis's old NR
@MODIFIED   :    version of mincfft
@MODIFIED   :
@MODIFIED   : Committing in .
@MODIFIED   :
@MODIFIED   : Added Files:
@MODIFIED   :  Makefile fft_support.c fft_support.h mincfft.c
@MODIFIED   : ----------------------------------------------------------------------
@MODIFIED   :
 * Revision 1.1  93/11/05  14:21:07  louis
 * Initial revision
 * 
---------------------------------------------------------------------------- */

#ifndef FFT_SUPPORT_H
#define FFT_SUPPORT_H


#include <volume_io.h>

#define   MAX_OUTFILES           8
#define   OUTPUT_REAL_AND_IMAG   0
#define   OUTPUT_REAL            1
#define   OUTPUT_IMAG            2
#define   OUTPUT_MAGNITUDE       3
#define   OUTPUT_MAGLN           4
#define   OUTPUT_MAG10           5
#define   OUTPUT_PHASE           6
#define   OUTPUT_POWER           7


Status prep_volume(Volume *tmp, Volume *data);
Status proj_volume(Volume *tmp, Volume *data, int job);
Status fft_volume(Volume data, int inverse_flg, int dim, int centre);



#endif
