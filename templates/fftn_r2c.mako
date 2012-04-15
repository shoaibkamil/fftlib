// fftn_R2C ... find FFT of in input real numpy array (time domain 
//              real or TR)returning a complex numpy array (frequency
//              domain complex or FC)
//

PyObject* SJSfftn_R2C(PyObject* TR) 
{

// Setup pointers ...
    double *dTR; //to data in the input numpy array. 
    fftw_complex *dFC; //to the output from FFTW

    dTR = (double*) PyArray_DATA(PyArray_GETCONTIGUOUS((PyArrayObject*)TR));

// I believe PyArray_DIMS will give me an array with the size 
// of each dimension. But I can't find anywhere a macro to tell 
// me how many dimensions there are!!!!  So I will import this
// as a template parameter (num_dims) for now. But I am sure 
// there is a way to extract this information from the Python Array.
//
    int i, ntot, nd;
    nd=${num_dims};
    int *ndimTR = (int*) malloc(nd*sizeof(int));
    npy_intp *ndimTR_np = (npy_intp*) PyArray_DIMS((PyArrayObject*)TR);
	for(i=0;i<nd;i++) ndimTR[i]=(int)ndimTR_np[i];
    for(i=0,ntot=1;i<nd;i++)ntot *= (*(ndimTR+i));

//
//  I need to look up the specific format required by the
//  fftw r2c.  I am setting aside too much space since 
//  symetry in the real to complex FFT can pack the
//  result into space for the real data (plus 1 .. I think)
//
// allocate space for the output, complex array from FFTW
// note: fftw defines complex as typedef double fftw_complex[2]
//       with element [0] real and [1] imaginary.
    
   dFC = (fftw_complex *)fftw_malloc(ntot * sizeof(fftw_complex)); 

// setup FFTW Plan (use r2c for arbitrary number of dimensions)
   fftw_plan pl_forward = fftw_plan_dft_r2c(nd, ndimTR, dTR,  
                                            dFC, FFTW_ESTIMATE);

// compute the fft
   fftw_execute(pl_forward);

//  
// I suspect I will need to create a new python object with the
// output data from FFTW which is the returned pyobject
//
// pack data into a numpy array object
//
    PyArrayObject* PyTmp = (PyArrayObject*)PyArray_SimpleNewFromData(nd,
                               ndimTR_np, NPY_CDOUBLE, dFC);
    Py_INCREF(PyTmp);
    
// free space and cleanup plans (TGM: I need to find a way to save 
// and reuse plans).   I need to think about the arrays I've created
// here and make sure I clean them up as well.  I am nervous about 
// doing so, however, since the returned python object uses the
// memory I've created for the FFT.

   fftw_destroy_plan(pl_forward);
   fftw_cleanup();

   return (PyObject*)PyTmp;
}
