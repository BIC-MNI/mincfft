/* fft_support.h */

#ifndef FFT_SUPPORT_H
#define FFT_SUPPORT_H


#include <minc2.h>
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


VIO_Status prep_volume(VIO_Volume *tmp, VIO_Volume *data);
VIO_Status proj_volume(VIO_Volume *tmp, VIO_Volume *data, int job);
VIO_Status fft_volume(VIO_Volume data, int inverse_flg, int dim, int centre);



#endif
