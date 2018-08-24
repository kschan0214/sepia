/*  kernel.h

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

// General kernel interpolation class

#if !defined(kernel_h)
#define kernel_h

#include <iostream>
#include <string>
#include <set>
#include <cmath>
#include "newmat.h"
      
using namespace NEWMAT;
using namespace std;

namespace MISCMATHS {

  /////////////////////////////////////////////////////////////////////////

  // Interpolation kernel storage class

  class kernelstorage
    {
    private:
      // NB: all widths are kernel half-widths (i.e. x \in [ -w, +w ] )
      int p_widthx;
      int p_widthy;
      int p_widthz;
      ColumnVector p_kernelx;
      ColumnVector p_kernely;
      ColumnVector p_kernelz;

      // stop all forms of creation except the constructors below
      kernelstorage();
      const kernelstorage& operator=(kernelstorage&);
      kernelstorage(kernelstorage&);


    public:
      float *storex;
      float *storey;
      float *storez;

      kernelstorage(const ColumnVector& kx, const ColumnVector& ky,
		    const ColumnVector& kz, int wx, int wy, int wz)	
        {
	  p_kernelx = kx; p_kernely = ky; p_kernelz = kz;
	  p_widthx = wx;  p_widthy = wy;  p_widthz = wz;
	  storez = new float[2*wz+1];
	  storey = new float[2*wy+1];
	  storex = new float[2*wx+1];
 	}

      ~kernelstorage()
      { 
	delete storex;
	delete storey;
	delete storez;
      }
      
      class comparer
	{
	public:
	  bool operator()(const kernelstorage* k1, 
			  const kernelstorage* k2) const
	    {
	      // comparison of sizes and values (toleranced)
	      if ( (k1->p_widthx!=k2->p_widthx) || 
		   (k1->p_widthy!=k2->p_widthy) || 
		   (k1->p_widthz!=k2->p_widthz) )
		return false;
	      if ( ( (k1->p_kernelx - k2->p_kernelx).MaximumAbsoluteValue()
		     > 1e-8 * k1->p_kernelx.MaximumAbsoluteValue() ) ||
		   ( (k1->p_kernely - k2->p_kernely).MaximumAbsoluteValue() 
		     > 1e-8 * k1->p_kernely.MaximumAbsoluteValue() ) ||
		   ( (k1->p_kernelz - k2->p_kernelz).MaximumAbsoluteValue() 
		     > 1e-8 * k1->p_kernelz.MaximumAbsoluteValue() ) )
		return false;
	      return true;
	    }
	};

      friend class comparer;

      int widthx() const { return p_widthx; }
      int widthy() const { return p_widthy; }
      int widthz() const { return p_widthz; }
      const ColumnVector& kernelx() const { return p_kernelx; }
      const ColumnVector& kernely() const { return p_kernely; }
      const ColumnVector& kernelz() const { return p_kernelz; }

    };


  /////////////////////////////////////////////////////////////////////////////

  class kernel
    {
    private:
      static set<kernelstorage*, kernelstorage::comparer> existingkernels;
      kernelstorage* storedkernel;

    public:
      kernel() { storedkernel = 0; }

      const kernel& operator=(const kernel& source)
      {
	// am allowed to copy pointers if other class either
	//  always exists or manages reference counts and self-deletes
	this->existingkernels = source.existingkernels;
	this->storedkernel = source.storedkernel;
	// signal storedkernel has an extra reference
	//   and that old storedkernel has one less reference
	return *this;
      }

      kernel(const kernel& source)
      {
	this->operator=(source);
      }

      virtual ~kernel() 
      { 
	// signal storedkernel it has one less reference
      }
      
      
      void setkernel(const ColumnVector& kx, const ColumnVector& ky,
		     const ColumnVector& kz, int wx, int wy, int wz)
      {		  
	// see if already in list:
	storedkernel = new kernelstorage(kx,ky,kz,wx,wy,wz);
	set<kernelstorage*, kernelstorage::comparer>::iterator 
	  it = existingkernels.find(storedkernel);
	if (it==existingkernels.end()) {		  
	  existingkernels.insert(storedkernel);
	  // signal that this is the first reference for storedkernel
	} else {
	  delete storedkernel;
	  storedkernel = *it;
	  // signal that *it has another reference now
	}
      }

      const kernelstorage* kernelvals() { return storedkernel; }
      
  };


  /////////////////////////////////////////////////////////////////////////

  //////// Support functions /////////
  
  float kernelval(float x, int w, const ColumnVector& kernel);
  float sincfn(float x);
  float hanning(float x, int w);
  float blackman(float x, int w);
  float rectangular(float x, int w);
  ColumnVector sinckernel1D(const string& sincwindowtype, int w, int n);
  kernel sinckernel(const string& sincwindowtype, int w, int nstore);
  kernel sinckernel(const string& sincwindowtype, 
		    int wx, int wy, int wz, int nstore);
  float extrapolate_1d(const ColumnVector& data, const int index);
  float interpolate_1d(const ColumnVector& data, const float index);
  float kernelinterpolation_1d(const ColumnVector& data, float index, const ColumnVector& userkernel, int width);
  float kernelinterpolation_1d(const ColumnVector& data, float index);
  float kernelinterpolation_1d(RowVector data, float index);
  float hermiteinterpolation_1d(const ColumnVector& data, int p1, int p4, float t);
}

#endif

