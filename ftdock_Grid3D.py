## INTEL CONFIDENTIAL
## Copyright 2005-2011 Intel Corporation.  All Rights Reserved.
##
## The source code contained or described herein and all documents related
## to the source code ("Material") are owned by Intel Corporation or its
## suppliers or licensors.  Title to the Material remains with Intel
## Corporation or its suppliers and licensors.  The Material contains
## trade secrets and proprietary and confidential information of Intel
## or its suppliers and licensors.  The Material is protected by worldwide
## copyright and trade secret laws and treaty provisions.  No part of the
## Material may be used, copied, reproduced, modified, published, uploaded,
## posted, transmitted, distributed, or disclosed in any way without
## Intel's prior express written permission.
##
## No license under any patent, copyright, trade secret or other
## intellectual property right is granted to or conferred upon you by
## disclosure or delivery of the Materials, either expressly, by
## implication, inducement, estoppel or otherwise.  Any license under such
## intellectual property rights must be express and approved by Intel in
## writing.

import numpy as np
import pyfftlib as fft

class Grid3D(object):
    '''Create 3D grid to discretize and correlate molecules'''

    def __init__(self, grid_dimension):
        '''Construct a 3D grid initialized to zero'''
        x = y = z = grid_dimension
        self.grid = np.zeros([x, y, z], dtype=float)

    def Transform(self):
        '''Return the FFT of the 3D grid'''
        return fft.lib().fftn_R2C(self.grid)
#       return np.fft.fftn(self.grid)

    def InverseTransform(self):
        '''Return the inverse FFT of the 3D grid'''
        return np.real_if_close(fft.lib().ifftn_C2R(self.grid), tol = 10000000)
#       return np.real_if_close(np.fft.ifftn(self.grid), tol = 10000000)

    def FFTShift(self):
        '''Return grid with half-spaces swapped for all axes'''
        return np.fft.fftshift(self.grid)

    def GetBestScore(self, idx):
        '''Return highest score in the correlation function'''
        return self.grid[idx]

    def GetBestTranslation(self):
        '''Return the highest scoring translation in the correlation function'''
        flat_idx = self.grid.argmax()
        dims = self.grid.shape
        return np.unravel_index(flat_idx, dims)   # hag: Debug, may need to convert to translation vector

##    def GetNBestTranslations(self, n):
##        print self.grid.argsort()

    def DiscretizeMolecule(self, molecule, interior_score, surface_score, surface_thickness, grid_resolution):
        '''Discretize a molecule in a 3D grid'''

        # Determine which grid nodes represent atoms
        for residue in molecule.structure.get_residues():
            full_id = residue.get_full_id()
            if (full_id[3][0] == ' '):   # Check hetero-flag, exclude heteroatoms
                for atom in residue:
                    (x, y, z) = atom.get_coord()
                    
                    # Convert atomic coordinates to upper and lower grid coordinates
                    i1 = int(round((x - 1.8) / grid_resolution))
                    i2 = int(round((x + 1.8) / grid_resolution)) + 1
                    j1 = int(round((y - 1.8) / grid_resolution))
                    j2 = int(round((y + 1.8) / grid_resolution)) + 1
                    k1 = int(round((z - 1.8) / grid_resolution))
                    k2 = int(round((z + 1.8) / grid_resolution)) + 1
                    self.grid[i1:i2,j1:j2,k1:k2] = surface_score

        # Find interior nodes of stationary molecule
        if (interior_score != surface_score):
            thickness = int(round(surface_thickness / grid_resolution))
            for (i,j,k) in np.transpose(self.grid.nonzero()):
                if np.all(self.grid[i-thickness:i+thickness+1,
                                    j-thickness:j+thickness+1,
                                    k-thickness:k+thickness+1] != 0.0):
                    self.grid[i,j,k] = interior_score

            # Verify that stationary molecule has some interior nodes
            if np.any(self.grid == interior_score):
                pass
            else:
                sys.exit('ftdock_Grid3D: warning: stationary molecule has no interior nodes, grid resolution probably too low')
