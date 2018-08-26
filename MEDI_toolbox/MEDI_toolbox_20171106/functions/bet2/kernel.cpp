/*  kernel.cc

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

#include "kernel.h"
#include "miscmaths.h"

namespace MISCMATHS {

  set<kernelstorage*, kernelstorage::comparer> kernel::existingkernels;  

  //////// Support function /////////

  float kernelval(float x, int w, const ColumnVector& kernel)
  {
    // linearly interpolates to get the kernel at the point (x)
    //   given the half-width w
    if (fabs(x)>w) return 0.0;
    float halfnk = (kernel.Nrows()-1.0)/2.0;
    float dn = x/w*halfnk + halfnk + 1.0;
    int n = (int) floor(dn);
    dn -= n;
    if (n>(kernel.Nrows()-1)) return 0.0;
    if (n<1) return 0.0;
    
    return kernel(n)*(1.0-dn) + kernel(n+1)*dn;
  }
  
  inline bool in_bounds(const ColumnVector& data, int index)
  { return ( (index>=1) && (index<=data.Nrows())); }
  
  inline bool in_bounds(const ColumnVector& data, float index)
  { return ( ((int)floor(index)>=1) && ((int)ceil(index)<=data.Nrows())); }
  
  float sincfn(float x)
  {
    if (fabs(x)<1e-7) { return 1.0-fabs(x); }
    float y=M_PI*x;
    return sin(y)/y;
  }
  
  float hanning(float x, int w)
  {  // w is half-width
    if (fabs(x)>w) 
      return 0.0;
    else
      return (0.5 + 0.5 *cos(M_PI*x/w));
  }
  
  float blackman(float x, int w)
  {  // w is half-width
    if (fabs(x)>w) 
      return 0.0;
    else
      return (0.42 + 0.5 *cos(M_PI*x/w) + 0.08*cos(2.0*M_PI*x/w));
  }
  
  float rectangular(float x, int w)
  {  // w is half-width
    if (fabs(x)>w) 
      return 0.0;
    else
      return 1.0;
  }

  ColumnVector sinckernel1D(const string& sincwindowtype, int w, int n)
  {  // w is full-width
    int nstore = n;
    if (nstore<1) nstore=1;
    ColumnVector ker(nstore);
    int hw = (w-1)/2; // convert to half-width
    // set x between +/- width
    float halfnk = (nstore-1.0)/2.0;
    for (int n=1; n<=nstore; n++) {
      float x=(n-halfnk-1)/halfnk*hw;
      if ( (sincwindowtype=="hanning") || (sincwindowtype=="h") ) {
	ker(n) = sincfn(x)*hanning(x,hw);
      } else if ( (sincwindowtype=="blackman") || (sincwindowtype=="b") ) {
	ker(n) = sincfn(x)*blackman(x,hw);
      } else if ( (sincwindowtype=="rectangular") || (sincwindowtype=="r") ) {
	ker(n) = sincfn(x)*rectangular(x,hw);
      } else {
	cerr << "ERROR: Unrecognised sinc window type - using Blackman" << endl;
	ker = sinckernel1D("b",w,nstore);
	return ker;
      }
    }
    return ker;
  }

  
  kernel sinckernel(const string& sincwindowtype, int w, int nstore) 
  {
    kernel sinck;
    sinck = sinckernel(sincwindowtype,w,w,w,nstore);
    return sinck;
  }


  kernel sinckernel(const string& sincwindowtype,
		    int wx, int wy, int wz, int nstore)
  { // widths are full-widths
    kernel sinckern;
    if (nstore<1) nstore=1;

    // convert all widths to half-widths
    int hwx = (wx-1)/2;
    int hwy = (wy-1)/2;
    int hwz = (wz-1)/2;

    ColumnVector kx, ky, kz;
    // calculate kernels
    kx = sinckernel1D(sincwindowtype,wx,nstore);
    ky = sinckernel1D(sincwindowtype,wy,nstore);
    kz = sinckernel1D(sincwindowtype,wz,nstore);
    
    sinckern.setkernel(kx,ky,kz,hwx,hwy,hwz);
    return sinckern;
  }

  // dummy fn for now
  float extrapolate_1d(const ColumnVector& data, const int index)
  {
    float extrapval;

    if (in_bounds(data, index))
      extrapval = data(index);
    else if (in_bounds(data, index-1))
      extrapval = data(data.Nrows());
    else if (in_bounds(data, index+1))
      extrapval = data(1);
    else
      extrapval = mean(data).AsScalar();

    return extrapval;
  }
  
  // basic trilinear call
  float interpolate_1d(const ColumnVector& data, const float index)
  {
    float interpval;
    int low_bound = (int)floor(index);
    int high_bound = (int)ceil(index); 

    if (in_bounds(data, index))
      interpval = data(low_bound) + (index - low_bound)*(data(high_bound) - data(low_bound));
    else
      interpval = extrapolate_1d(data, round(index));

    return interpval;
  }

    
  //////// Spline Support /////////

  float hermiteinterpolation_1d(const ColumnVector& data, int p1, int p4, float t)
  {
    // Q(t) = (2t^3 - 3t^2 + 1)P_1 + (-2t^3 + 3t^2)P_4 + (t^3 - 2t^2 + t)R_1 + (t^3 - t^2)R_4
    // inputs: points P_1, P_4; tangents R_1, R_4; interpolation index t;

    float retval, r1 = 0.0, r4 = 0.0;
    
    if (!in_bounds(data,p1) || !in_bounds(data,p4)) {
      cerr << "Hermite Interpolation - ERROR: One or more indicies lie outside the data range. Returning ZERO" << endl;
      retval = 0.0;
    } else if ((t < 0) || (t > 1)) {
      cerr << "Hermite Interpolation - ERROR: Interpolation index must lie between 0 and 1. Returning ZERO" << endl;
      retval = 0.0;
      /*    } else if (t == 0.0) {
	    retval = data(p1);
	    } else if (t == 1.0) {
	    retval = data(p4);
      */   
    } else {
      r1 = 0.5 * (extrapolate_1d(data, p1) - extrapolate_1d(data, p1 - 1)) + 0.5 * (extrapolate_1d(data, p1 + 1) - extrapolate_1d(data, p1));// tangent @ P_1
      r4 = 0.5 * (extrapolate_1d(data, p4) - extrapolate_1d(data, p4 - 1)) + 0.5 * (extrapolate_1d(data, p4 + 1) - extrapolate_1d(data, p4));// tangent @ P_4
      
      float t2 = t*t; float t3 = t2*t;
      retval = (2*t3 - 3*t2 + 1)*data(p1) + (-2*t3 + 3*t2)*data(p4) + (t3 - 2*t2 + t)*r1 + (t3 - t2)*r4;
    }

    // cerr << "p1, p4, t, r1, r4 = " << p1 << ", " << p4 << ", " << t << ", " << r1 << ", " << r4 << endl;
    
    return retval;
  }


  //////// Kernel Interpolation Call /////////

  float kernelinterpolation_1d(const ColumnVector& data, float index, const ColumnVector& userkernel, int width)
  {
    int widthx = (width - 1)/2;
    // kernel half-width  (i.e. range is +/- w)
    int ix0;
    ix0 = (int) floor(index);
    
    int wx(widthx);
    vector<float> storex(2*wx+1);
    for (int d=-wx; d<=wx; d++) 
       storex[d+wx] = kernelval((index-ix0+d),wx,userkernel);

    float convsum=0.0, interpval=0.0, kersum=0.0;
    
    int xj;
    for (int x1=ix0-wx; x1<=ix0+wx; x1++) {
      if (in_bounds(data, x1)) {
	xj=ix0-x1+wx;
	float kerfac = storex[xj];
	convsum += data(x1) * kerfac;
	kersum += kerfac;
      }
    }

    if ( (fabs(kersum)>1e-9) ) {
      interpval = convsum / kersum;
    } else {
      interpval = (float) extrapolate_1d(data, ix0);
    }
    return interpval;
  }
  
  ////// Kernel wrappers //////

  float kernelinterpolation_1d(const ColumnVector& data, float index)
  {
    ColumnVector userkernel = sinckernel1D("hanning", 7, 1201);
    return kernelinterpolation_1d(data, index, userkernel, 7);
  }

  float kernelinterpolation_1d(RowVector data, float index)
  {
    ColumnVector userkernel = sinckernel1D("hanning", 7, 1201);
    return kernelinterpolation_1d(data.t(), index, userkernel, 7);
  }
}

