#include <stdio.h>
#include <complex.h>
#include <fftw3.h>

int main()
{
   int i, grd_size = 16;
   double *fA, *fB, *fC;
   fftw_complex *FA, *FB, *FC;
   fftw_plan fA_to_FA, fB_to_FB, FC_to_fC;

   /* Allocate grids representing discretized molecules */
   fA = (double *) fftw_malloc (sizeof(double) * grd_size);
   fB = (double *) fftw_malloc (sizeof(double) * grd_size);
   fC = (double *) fftw_malloc (sizeof(double) * grd_size);

   FA = (fftw_complex *) fftw_malloc (sizeof(fftw_complex) * grd_size);
   FB = (fftw_complex *) fftw_malloc (sizeof(fftw_complex) * grd_size);
   FC = (fftw_complex *) fftw_malloc (sizeof(fftw_complex) * grd_size);

   /* Initialize grids representing discretized molecules */
   for (i = 0; i < grd_size; i++)
   {
      fA[i] = fB[i] = 0.0;
   }
   fA[6] = fB[6] = 1.0;
   fA[7] = fB[7] = 1.0;
   fA[8] = fB[8] = 1.0;
   fA[9] = fB[9] = 1.0;
   fA[7] = -5.0;

   /* Transform "molecule" grids, note that FFTW computes unnormalized transforms */
   fA_to_FA = fftw_plan_dft_r2c_1d(grd_size, fA, FA, FFTW_MEASURE);
   fB_to_FB = fftw_plan_dft_r2c_1d(grd_size, fB, FB, FFTW_MEASURE);
   fftw_execute(fA_to_FA);
   fftw_execute(fB_to_FB);

   printf("FA = ");
   for (i = 0; i < grd_size; i++) printf ("(%.1f %.1fi) ", creal(FA[i]), cimag(FA[i]));
   printf("\n\n");
   printf("FB = ");
   for (i = 0; i < grd_size; i++) printf ("(%.1f %.1fi) ", creal(FB[i]), cimag(FB[i]));
   printf("\n\n");

   /* Perform Fourier correlation */
   for (i = 0; i < grd_size; i++)
   {
      FC[i] = conj(FA[i]) * FB[i];
   }

   /* Perform inverse transform */
   FC_to_fC = fftw_plan_dft_c2r_1d(grd_size, FC, fC, FFTW_MEASURE);
   fftw_execute(FC_to_fC);

   /* Normalize the correlation function */
   for (i = 0; i < grd_size; i++) fC[i] = fC[i] / grd_size;

   printf("fC = ");
   for (i = 0; i < grd_size; i++) printf ("%.1f ", fC[i]);
   printf("\n");

   /* Cleanup */
   fftw_destroy_plan(fA_to_FA);
   fftw_destroy_plan(fB_to_FB);
   fftw_destroy_plan(FC_to_fC);
   fftw_free(fA); fftw_free(fB); fftw_free(fC);
   fftw_free(FA); fftw_free(FB); fftw_free(FC);
}
