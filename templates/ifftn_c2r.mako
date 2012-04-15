// ifftn_C2R ... find inverse FFT of in input complex numpy array 
//              (frequency domain complex or FC) and return a real 
//              numpy array time domain real or TR).
//
//              Note: the complex array FC  is assumed to have been
//              produced by finding the forward transform of a
//              real time domain array.


PyObject* SJSifftn_C2R(PyObject* FC) 
{

// Setup pointers ...
    fftw_complex *dFC; //to data in the input numpy array. 
    double *dTR; //to the output from FFTW

    dFC = (fftw_complex*) PyArray_DATA(PyArray_GETCONTIGUOUS((PyArrayObject*)FC));

// I believe PyArray_DIMS will give me an array with the size 
// of each dimension. But I can't find anywhere a macro to tell 
// me how many dimensions there are!!!!  So I will import this
// as a template parameter (num_dims) for now. But I am sure 
// there is a way to extract this information from the Python Array.
//
    int i, ntot, nd;
    nd=${num_dims};
    int *ndimFC = (int*) malloc(nd*sizeof(int));
    npy_intp *ndimFC_np = (npy_intp*) PyArray_DIMS((PyArrayObject*)FC);
    for(i=0;i<nd;i++) ndimFC[i]=(int)ndimFC_np[i];
    for(i=0,ntot=1;i<nd;i++)ntot *= (*(ndimFC+i));


//
// allocate space for the output, double array from FFTW
    
   dTR = (double *)fftw_malloc(ntot * sizeof(double)); 

// setup FFTW Plan (use c2r for arbitrary number of dimensions)
   fftw_plan pl_inverse = fftw_plan_dft_c2r(nd, ndimFC, dFC,  
                                            dTR, FFTW_ESTIMATE);

// compute the fft
   fftw_execute(pl_inverse);

//  
// I suspect I will need to create a new python object with the
// output data from FFTW which is the returned pyobject
//
// pack data into a numpy array object
//
    PyArrayObject* PyTmp = (PyArrayObject*)PyArray_SimpleNewFromData(nd,
                               ndimFC_np, NPY_DOUBLE, dTR);
    Py_INCREF(PyTmp);
    
// free space and cleanup plans (TGM: I need to find a way to save 
// and reuse plans).   I need to think about the arrays I've created
// here and make sure I clean them up as well.  I am nervous about 
// doing so, however, since the returned python object uses the
// memory I've created for the FFT.

   fftw_destroy_plan(pl_inverse);
   fftw_cleanup();

   return (PyObject*)PyTmp;
}

