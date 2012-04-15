import numpy as np
import pyfftlib as fft
np.set_printoptions(precision = 1, linewidth = 80)
grd_size = 16

# Initialize 1D grids
fA = np.zeros([grd_size], dtype=float)
FA = np.zeros([grd_size], dtype=complex)
fB = np.zeros([grd_size], dtype=float)
FB = np.zeros([grd_size], dtype=complex)
fC = np.zeros([grd_size], dtype=float)
FC = np.zeros([grd_size], dtype=complex)

# Setup 1D grids representing discretized molecules
# fA = (0,  0,  0,  0,  0,  0,  1, -5,  1,  1,  0,  0,  0,  0,  0,  0)
# fB = (0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  0,  0,  0,  0,  0,  0)
fA[6:10] = 1.0
fA[7] = -5.0
fB[6:10] = 1.0
print
print 'fA =', fA.round()
print 'fB =', fB.round()
print

# Transform "molecule" grids
FA = fft.lib().fftn_R2C(fA)
FB = fft.lib().fftn_R2C(fB)
#FA = np.fft.fft(fA)
#FB = np.fft.fft(fB)

# Perform Fourier correlation
FC = FA.conj() * FB

# Perform inverse transform
fC = np.real_if_close(fft.lib().ifftn_C2R(FC))

# Look at the correlation function (fC) and the position of the best score
# before shifting.
print 'fC =', fC.round(), 'no shift'
print 'fC max =', fC.max(), 'fC maxloc =', fC.argmax()
print

# It is not obvious how the best correlation score in fC relates to the
# translation vector for molecule B (fB) relative to molecule A (fA) so a
# half-interval shift is performed so that the optimum translation vector
# is determined relative to grid center, which is the origin in Cartesian
# space.

fC = np.fft.fftshift(fC)
print 'fC =', fC.round(), 'half-interval shift'
print 'fC max =', fC.max(), 'fC maxloc =', fC.argmax()
print

# The output is now more intuitive for molecular docking. The negative
# correlation scores (i.e., the molecules overlap too much) are in the
# center of the grid, the zeroes (i.e., the molecules are not in contact)
# are on the edges, and the positive score (i.e., the molecules are in
# contact without excessive penetration) are in between. In this example,
# the translation vector that gives the optimum correlation score is simply
# [grd_size/2 - fC.argmax()] = 2 so fB is shifted two grid nodes to the right:
#
# fA = (0,  0,  0,  0,  0,  0,  1, -5,  1,  1,  0,  0,  0,  0,  0,  0)
# fB =         (0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  0,  0,  0,  0,  0,  0)
#
# In Cartesian space, (2 * grid_resolution) equals the number of Angstroms
# to tranlate molecule B on the x-axis.
