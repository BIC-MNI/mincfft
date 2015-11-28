/* mincfft.c                                                                 */
/*                                                                           */
/* A widget for FFT'ing MINC data                                            */
/*                                                                           */
/* Louis Collins- louis@bic.mni.mcgill.ca                                    */
/* Brain Imaging Centre                                                      */
/* Montreal Neurological Institute                                           */
/*                                                                           */
/* Andrew Janke - a.janke@gmail.com                                          */
/* Center for Magnetic Resonance                                             */
/* University of Queensland                                                  */
/*                                                                           */
/* Copyright Andrew Janke, The University of Queensland & Louis Collins,     */
/* Montreal Neurological Institute, Canada.                                  */
/* Permission to use, copy, modify, and distribute this software and its     */
/* documentation for any purpose and without fee is hereby granted,          */
/* provided that the above copyright notice appear in all copies.  The       */
/* author and the University of Queensland make no representations about the */
/* suitability of this software for any purpose.  It is provided "as is"     */
/* without express or implied warranty.                                      */


#include <float.h>
#include <ParseArgv.h>
#include <time_stamp.h>
#include "fft_support.h"

/* function prototypes */
void print_version_info(void);

/* hack for pretty-printing */
char *out_names[MAX_OUTFILES] = {
   "real+imag",
   "real     ",
   "imag     ",
   "magnitude",
   "magln    ",
   "mag10    ",
   "phase    ",
   "power    "
};

int      verbose = FALSE;
int      clobber = FALSE;
int      inv_fft = FALSE;
int      centre_fft = FALSE;
int      fft_dim = 3;
char    *outfiles[MAX_OUTFILES] = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL };
int      is_signed = FALSE;
nc_type  dtype = NC_FLOAT;

static ArgvInfo argTable[] = {
   {NULL, ARGV_HELP, (char *)NULL, (char *)NULL,
    "General options:"},
   {"-version", ARGV_FUNC, (char *)print_version_info, (char *)NULL,
    "print version info and exit"},
   {"-verbose", ARGV_CONSTANT, (char *)TRUE, (char *)&verbose,
    "Print out extra information."},
   {"-clobber", ARGV_CONSTANT, (char *)TRUE, (char *)&clobber,
    "Clobber existing files."},

   {NULL, ARGV_HELP, NULL, NULL, "\nOutfile Options"},
   {"-byte", ARGV_CONSTANT, (char *)NC_BYTE, (char *)&dtype,
    "Write out byte data."},
   {"-short", ARGV_CONSTANT, (char *)NC_SHORT, (char *)&dtype,
    "Write out short integer data."},
   {"-long", ARGV_CONSTANT, (char *)NC_LONG, (char *)&dtype,
    "Write out long integer data."},
   {"-float", ARGV_CONSTANT, (char *)NC_FLOAT, (char *)&dtype,
    "Write out single-precision data. (Default)"},
   {"-double", ARGV_CONSTANT, (char *)NC_DOUBLE, (char *)&dtype,
    "Write out double-precision data."},
   {"-signed", ARGV_CONSTANT, (char *)TRUE, (char *)&is_signed,
    "Write signed integer data."},
   {"-unsigned", ARGV_CONSTANT, (char *)FALSE, (char *)&is_signed,
    "Write unsigned integer data."},

   {NULL, ARGV_HELP, NULL, NULL, "\nFFT options"},
   {"-2D", ARGV_CONSTANT, (char *)2, (char *)&fft_dim,
    "Do a 2D FFT (Default 3D)."},
   {"-forward", ARGV_CONSTANT, (char *)FALSE, (char *)&inv_fft,
    "Calculate the forward FFT (default)."},
   {"-inverse", ARGV_CONSTANT, (char *)TRUE, (char *)&inv_fft,
    "Calculate the inverse FFT."},
   {"-centre", ARGV_CONSTANT, (char *)TRUE, (char *)&centre_fft,
    "Re-orient quadrants to force resulting data to the centre"},
   {"-center", ARGV_CONSTANT, (char *)TRUE, (char *)&centre_fft,
    "Synonym for our North American friends"},

   {NULL, ARGV_HELP, NULL, NULL, "\nOutput file types for FFT"},
   {"-both", ARGV_STRING, (char *)1, (char *)&outfiles[OUTPUT_REAL_AND_IMAG],
    "<file.mnc> Complex Real and Imaginary data (default)."},
   {"-real", ARGV_STRING, (char *)1, (char *)&outfiles[OUTPUT_REAL],
    "<file.mnc> Real component of data."},
   {"-imaginary", ARGV_STRING, (char *)1, (char *)&outfiles[OUTPUT_IMAG],
    "<file.mnc> Imaginary component of data."},
   {"-magnitude", ARGV_STRING, (char *)1, (char *)&outfiles[OUTPUT_MAGNITUDE],
    "<file.mnc> magnitude of Real and Imaginary data."},
   {"-magln", ARGV_STRING, (char *)1, (char *)&outfiles[OUTPUT_MAGLN],
    "<file.mnc> ln magnitude of Real and Imaginary data."},
   {"-mag10", ARGV_STRING, (char *)1, (char *)&outfiles[OUTPUT_MAG10],
    "<file.mnc> log 10 magnitude of Real and Imaginary data."},
   {"-phase", ARGV_STRING, (char *)1, (char *)&outfiles[OUTPUT_PHASE],
    "<file.mnc> phase of Real and Imaginary data."},
   {"-power", ARGV_STRING, (char *)1, (char *)&outfiles[OUTPUT_POWER],
    "<file.mnc> power spectrum."},

   {NULL, ARGV_HELP, NULL, NULL, ""},
   {NULL, ARGV_END, NULL, NULL, NULL}
};

char *spac_dimorder[] = { MIzspace, MIyspace, MIxspace };
char *freq_dimorder[] = { MIzspace, MIyspace, MIxspace, MIvector_dimension };

int main(int argc, char *argv[]){
   char *in_fn;
   char *history;
   VIO_Status status;
   VIO_Volume tmp;
   VIO_Volume data;
   VIO_Volume *vol_ptr = NULL;
   int      c, in_ndims;
   int      n_outfiles;
   VIO_Real     min, max;

   minc_input_options in_ops;

   /* get the history string */
   history = time_stamp(argc, argv);

   /* get args */
   if(ParseArgv(&argc, argv, argTable, 0) || (argc < 2)){
      fprintf(stderr,
              "\nUsage: %s [<options>] <infile.mnc> [-outtype <type.mnc>] [<outfile.mnc>]\n",
              argv[0]);
      fprintf(stderr, "       %s [-help]\n\n", argv[0]);
      exit(EXIT_FAILURE);
      }
   in_fn = argv[1];
   if(argc > 2){
      outfiles[OUTPUT_REAL_AND_IMAG] = argv[2];
      }

   /* check for infile and outfiles */
   if(!file_exists(in_fn)){
      fprintf(stderr, "%s: Couldn't find input file %s.\n", argv[0], in_fn);
      exit(EXIT_FAILURE);
      }
   n_outfiles = 0;
   for(c = 0; c < MAX_OUTFILES; c++){
      if(outfiles[c] != NULL){
         if(!clobber && file_exists(outfiles[c])){
            fprintf(stderr, "%s: File %s exists, use -clobber to overwrite.\n", argv[0],
                    outfiles[c]);
            exit(EXIT_FAILURE);
            }
         n_outfiles++;
         }
      }
   if(n_outfiles == 0){
      fprintf(stderr, "%s: You should specify at least one outfile!\n", argv[0]);
      exit(EXIT_FAILURE);
      }


   /* read in the input file */
   in_ndims = get_minc_file_n_dimensions(in_fn);
   set_default_minc_input_options(&in_ops);
   set_minc_input_vector_to_scalar_flag(&in_ops, FALSE);
   if(in_ndims == 4){
      status = input_volume(in_fn, 4, freq_dimorder,
                            NC_UNSPECIFIED, FALSE, 0.0, 0.0, TRUE, &data, &in_ops);
      }
   else{
      status = input_volume(in_fn, 3, spac_dimorder,
                            NC_UNSPECIFIED, FALSE, 0.0, 0.0, TRUE, &tmp, &in_ops);
      status &= prep_volume(&tmp, &data);
      delete_volume(tmp);
      }

   if(status != VIO_OK){
      fprintf(stderr, "Problems reading: %s\n", in_fn);
      exit(EXIT_FAILURE);
      }

   if(verbose){
      VIO_Real     min_value, max_value;

      get_volume_real_range(data, &min_value, &max_value);

      fprintf(stdout, " | Input file:     %s\n", in_fn);
      fprintf(stdout, " | Input ndims:    %d\n", in_ndims);
      fprintf(stdout, " | min/max:        [%8.3f:%8.3f]\n", min_value, max_value);
      fprintf(stdout, " | Output files:\n");
      for(c = 0; c < MAX_OUTFILES; c++){
         if(outfiles[c] != NULL){
            fprintf(stdout, " |   [%d]:         %s => %s\n", c, out_names[c],
                    outfiles[c]);
            }
         }
      fprintf(stdout, " | FFT order:      %d\n", fft_dim);
      }

   /* FFT the volume */
   if(fft_volume(data, inv_fft, fft_dim, centre_fft) != VIO_OK){
      print_error("Problems during FFT of: %s", in_fn);
      }

   /* output the resulting volume(s) */
   for(c = 0; c < MAX_OUTFILES; c++){

      if(outfiles[c] != NULL){
         if(verbose){
            fprintf(stdout, "Outputting %s (%s) \t=> ", out_names[c], outfiles[c]);
            fflush(stdout);
            }

         /* do the projection if neccesarry */
         tmp = NULL;
         if(c == OUTPUT_REAL_AND_IMAG){
            VIO_Real     value;
            int      i, j, k, l;
            int      sizes[4];

            get_volume_sizes(data, sizes);

            /* set up max and min values */
            min = DBL_MAX;
            max = -DBL_MAX;
            for(i = sizes[0]; i--;){
               for(j = sizes[1]; j--;){
                  for(k = sizes[2]; k--;){
                     for(l = sizes[3]; l--;){

                        GET_VALUE_4D(value, data, i, j, k, l);
                        if(value > max){
                           max = value;
                           }
                        if(value < min){
                           min = value;
                           }
                        }
                     }
                  }
               }

            set_volume_real_range(data, min, max);

            vol_ptr = &data;
            }
         else{
            status = proj_volume(&data, &tmp, c);
            vol_ptr = &tmp;
            }

         if(verbose){
            get_volume_real_range(*vol_ptr, &min, &max);
            fprintf(stdout, "| range: [%g:%g]\n", min, max);
            }

         if(output_modified_volume(outfiles[c],
                                    dtype, is_signed, 0, 0,
                                    *vol_ptr, in_fn, history, NULL) != VIO_OK){
            print_error("Problems outputing: %s", outfiles[c]);
            }

         if(tmp != NULL){
            delete_volume(tmp);
            }
         }
      }


   delete_volume(data);
   return (status);
   }

void print_version_info(void){
   fprintf(stdout, "%s version %s\n", PACKAGE_NAME, PACKAGE_VERSION);
   fprintf(stdout, "Comments to %s\n", PACKAGE_BUGREPORT);
   fprintf(stdout, "\n");
   exit(EXIT_SUCCESS);
   }
