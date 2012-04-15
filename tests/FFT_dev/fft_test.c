//--------------------------------------------------------------
//
// program to test the FFT routines in FFTW
//
// I am testing the real to real routines which use a 
// half-complex format to hold forward transform results
// (see FFTW documentation for details).
//
// FFTW carries out un-normalized FFT's so the result from
// forward followed by an inverse transforms is scalled by N.
//
//--------------------------------------------------------------

#include <fftw3.h>
#include <math.h>
#define TwoPi  6.28318531
#define DEBUG  1    // non zero for debug output, 5 for lots of it
#define N      32   // length of FFTs 

//--------------------------------------------------------------
//  output an array printing only a few significant figures
//--------------------------------------------------------------
void out_arr(double *arr)
{
  for(int i=0;i<N;i++){
       if(4*(i/4) == i) printf("\n");  
       printf(" %10.3f ",*(arr+i));
  }
}

//--------------------------------------------------------------
//  Test the FFT ... make sure final scaled time domain is
//  similar to the orignial data. Returns the root mean 
//  square error.
//--------------------------------------------------------------
double testFFT(double *orig, double *final)
{
  double dif, errsq=0.0;
  double dN = (double)N;
  for(int i=0;i<N;i++){
       dif = (*(orig+i))*dN  - (*(final+i)) ;
 if(DEBUG == 5)printf("\n %f %f %f ", *(orig+i), *(final+i), dif); 
       errsq += dif*dif;
  }
  return sqrt(errsq);
}

//--------------------------------------------------------------
int main()
{
  double *tdom_in;   // orignial input time sequence
  double *fdom_out;  // half complex frequency domain
  double *tdom_out;  // return to time domain to recover tdom_in*N

  // parameters to produce time sequence with known spectrum
  double k = TwoPi/(double)N;
  double freq1 = 2.0, freq2 = 5.0,  freq3 = 8.0;
  double tcon = 0.25;
  double dN = (double)N;
  
  printf("test the fft software \n");
  
  // allocate data ... FFTW malloc produces alligned arrays
  tdom_in  = (double*) fftw_malloc(N*sizeof(double));
  fdom_out = (double*) fftw_malloc(N*sizeof(double));
  tdom_out = (double*) fftw_malloc(N*sizeof(double));

  // generate data with a known spectrum
  for(int i=0;i<N;i++){
     *(tdom_in+i) = tcon + (0.5*sin(freq1*(double)i*k) 
                          * 0.75*sin(freq2*(double)i*k)
                          * 0.25*sin(freq3*(double)i*k));
  }
  
  if(DEBUG){printf("\n orig data");  out_arr(tdom_in);}

  // setup plans for forward and inverse real to real FFTs where
  // we estimate rather than measure the best algorithm to use
  fftw_plan pl_forward = fftw_plan_r2r_1d(
                 N, tdom_in,  fdom_out, FFTW_R2HC, FFTW_ESTIMATE);

  fftw_plan pl_reverse = fftw_plan_r2r_1d(
                 N, fdom_out, tdom_out, FFTW_HC2R, FFTW_ESTIMATE);

  // Execute the FFTs
  fftw_execute(pl_forward);
  if(DEBUG){printf("\n freq domain");  out_arr(fdom_out);}

  fftw_execute(pl_reverse);
  if(DEBUG){printf("\n final time domain");  out_arr(tdom_out);}

  printf("\nerror = %f \n", testFFT(tdom_in, tdom_out));
  
  // free space and clean up plans
  fftw_destroy_plan(pl_forward);     
  fftw_destroy_plan(pl_reverse);
  fftw_cleanup();
  fftw_free(tdom_in);
  fftw_free(fdom_out);
  fftw_free(tdom_out);

  printf("the tests are done\n");
}
    
