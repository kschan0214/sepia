
import os
import glob

import nibabel as nib
import numpy as np
import scipy.io
from pathlib import Path

'''

prepare training dataset.

'''

# input_fn = Path('/testing/testing_lp_cnn/dipole.mat')
# output_fn = Path('/testing/testing_lp_cnn/dipole.npy')

dipole = scipy.io.loadmat(input_fn)['C']
dipole = np.swapaxes(dipole, 0, 1)

np.save(output_fn, dipole)
print('dipole saved.')

