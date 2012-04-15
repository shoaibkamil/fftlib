#
# The python version of fftlib.
#
# This code just calls the numpy fft routines.
#

import numpy as np

class lib(object):
    
    def __init__(self):
        self.pure_python = True

    def fftn_R2C(self, x_in):
        '''Return the N-dim DFT (complex) of a real array'''
        return np.fft.fftn(x_in)

    def ifftn_C2R(self,x_in):
        '''Return the N-dim inv DFT (real) of a complex array'''
        return np.real_if_close(np.fft.ifftn(x_in), tol = 10000000)
