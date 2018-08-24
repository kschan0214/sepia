/*  imfft.h

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 2001 University of Oxford  */

/*  Part of FSL - FMRIB's Software Library
    http://www.fmrib.ox.ac.uk/fsl
    fsl@fmrib.ox.ac.uk
    
    Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
    Imaging of the Brain), Department of Clinical Neurology, Oxford
    University, Oxford, UK
    
    
    LICENCE
    
    FMRIB Software Library, Release 4.0 (c) 2007, The University of
    Oxford (the "Software")
    
    The Software remains the property of the University of Oxford ("the
    University").
    
    The Software is distributed "AS IS" under this Licence solely for
    non-commercial use in the hope that it will be useful, but in order
    that the University as a charitable foundation protects its assets for
    the benefit of its educational and research purposes, the University
    makes clear that no condition is made or to be implied, nor is any
    warranty given or to be implied, as to the accuracy of the Software,
    or that it will be suitable for any particular purpose or for use
    under any specific conditions. Furthermore, the University disclaims
    all responsibility for the use which is made of the Software. It
    further disclaims any liability for the outcomes arising from using
    the Software.
    
    The Licensee agrees to indemnify the University and hold the
    University harmless from and against any and all claims, damages and
    liabilities asserted by third parties (including claims for
    negligence) which arise directly or indirectly from the use of the
    Software or the sale of any products based on the Software.
    
    No part of the Software may be reproduced, modified, transmitted or
    transferred in any form or by any means, electronic or mechanical,
    without the express permission of the University. The permission of
    the University is not required if the said reproduction, modification,
    transmission or transference is done without financial return, the
    conditions of this Licence are imposed upon the receiver of the
    product, and all original and amended source code is included in any
    transmitted product. You may be held legally responsible for any
    copyright infringement that is caused or encouraged by your failure to
    abide by these terms and conditions.
    
    You are not permitted under this Licence to use this Software
    commercially. Use for which any financial return is received shall be
    defined as commercial use, and includes (1) integration of all or part
    of the source code or the Software into a product for sale or license
    by or on behalf of Licensee to third parties or (2) use of the
    Software or any derivative of it for research with the final aim of
    developing software products for sale or license to a third party or
    (3) use of the Software or any derivative of it for research with the
    final aim of developing non-software products for sale or license to a
    third party, or (4) use of the Software to provide any service to an
    external organisation for which payment is received. If you are
    interested in using the Software commercially, please contact Isis
    Innovation Limited ("Isis"), the technology transfer company of the
    University, to negotiate a licence. Contact details are:
    innovation@isis.ox.ac.uk quoting reference DE/1112. */


// FFT routines for 3D images
#if !defined(__imfft_h)
#define __imfft_h

#include <string>
#include <iostream>
#include <fstream>
#if defined(_WIN32) || defined(_WIN64)
#include <io.h>
#else
#include <unistd.h>
#endif

#include "newmatap.h"
#include "newmatio.h"
#include "newimage.h"
#include "complexvolume.h"

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

using namespace MISCMATHS;
using namespace NEWMAT;

namespace NEWIMAGE {

////////////////////////////////////////////////////////////////////////////

int ifft2(complexvolume& vol);
int ifft2(const complexvolume& vin, complexvolume& vout);

int fft2(complexvolume& vol);
int fft2(const complexvolume& vin, complexvolume& vout);

int ifft3(complexvolume& vol);
int ifft3(const complexvolume& vin, complexvolume& vout);

int fft3(complexvolume& vol);
int fft3(const complexvolume& vin, complexvolume& vout);

void fftshift(complexvolume& vol);

//------------------------------------------------------------------------//

int ifft2(volume<float>& realvol, volume<float>& imagvol);
int ifft2(const volume<float>& realvin, const volume<float>& imagvin,
	  volume<float>& realvout, volume<float>& imagvout);

int fft2(volume<float>& realvol, volume<float>& imagvol);
int fft2(const volume<float>& realvin, const volume<float>& imagvin,
	 volume<float>& realvout, volume<float>& imagvout);

int ifft3(volume<float>& realvol, volume<float>& imagvol);
int ifft3(const volume<float>& realvin, const volume<float>& imagvin,
	  volume<float>& realvout, volume<float>& imagvout);

int fft3(volume<float>& realvol, volume<float>& imagvol);
int fft3(const volume<float>& realvin, const volume<float>& imagvin,
	 volume<float>& realvout, volume<float>& imagvout);

void fftshift(volume<float>& vol);

//------------------------------------------------------------------------//

int ifft2(volume<double>& realvol, volume<double>& imagvol);
int ifft2(const volume<double>& realvin, const volume<double>& imagvin,
	  volume<double>& realvout, volume<double>& imagvout);

int fft2(volume<double>& realvol, volume<double>& imagvol);
int fft2(const volume<double>& realvin, const volume<double>& imagvin,
	 volume<double>& realvout, volume<double>& imagvout);

int ifft3(volume<double>& realvol, volume<double>& imagvol);
int ifft3(const volume<double>& realvin, const volume<double>& imagvin,
	  volume<double>& realvout, volume<double>& imagvout);

int fft3(volume<double>& realvol, volume<double>& imagvol);
int fft3(const volume<double>& realvin, const volume<double>& imagvin,
	 volume<double>& realvout, volume<double>& imagvout);

void fftshift(volume<double>& vol);

////////////////////////////////////////////////////////////////////////////

}

#endif

