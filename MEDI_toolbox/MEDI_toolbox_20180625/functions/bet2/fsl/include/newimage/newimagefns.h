/*  newimagefns.h

    Mark Jenkinson et al, FMRIB Image Analysis Group

    Copyright (C) 2000-2006 University of Oxford  */

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

// General image processing functions

#if !defined(__newimagefns_h)
#define __newimagefns_h

#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include "newimage.h"
#include "newmatio.h"
#include "fslio/fslio.h"
#include "miscmaths/miscmaths.h"
#include "complexvolume.h"
#include "imfft.h"
#include <queue>

#ifndef MAX
#define MAX(a,b) (((a)>(b))?(a):(b))
#endif

#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif

using namespace std;
using namespace NEWMAT;

namespace NEWIMAGE {


  // The following lines are ignored by the current SGI compiler
  //  (version egcs-2.91.57)
  // A temporary fix of including the std:: in front of all abs() etc
  //  has been done for now
  using std::abs;
  using std::sqrt;

  ///////////////////////////////////////////////////////////////////////////

  // DIAGNOSTICS

  template <class T>
  void print_size(const volume<T>& source);

  template <class T>
  void print_volume_info(const volume<T>& source, const string& name="");
  template <class T>
  void print_volume_info(const volume4D<T>& source, const string& name="");

  ///////////////////////////////////////////////////////////////////////////

  // BASIC IMAGE SUPPORT FUNCTIONS

  // Complex volume support
  volume<float> abs(const volume<float>& realvol, 
		    const volume<float>& imagvol);

  volume<float> phase(const volume<float>& realvol, 
		      const volume<float>& imagvol);

  volume<float> real(const volume<float>& absvol, 
		     const volume<float>& phasevol);

  volume<float> imag(const volume<float>& absvol, 
		     const volume<float>& phasevol);

  // Basic Arithmetic Operations
  template <class T>
  volume<T> abs(const volume<T>& vol);
  
  template <class T>
  volume4D<T> abs(const volume4D<T>& vol);

  template <class T, class S>
  volume<T> divide(const volume<T>& numervol, const volume<T>& denomvol,
		   const volume<S>& mask);
  template <class T, class S>
  volume4D<T> divide(const volume4D<T>& numervol, const volume4D<T>& denomvol,
		     const volume<S>& mask);
  template <class T, class S>
  volume4D<T> divide(const volume4D<T>& numervol, const volume<T>& denomvol,
		     const volume<S>& mask);

  
  template <class T, class S>
  volume<T> mask_volume( const volume<T>& invol,const volume<S>& mask );

  template<class T>
  void indexadd(volume<T>& vola, const volume<T>& volb,
		const Matrix& indices);
  template<class T>
  void indexsubtract(volume<T>& vola, const volume<T>& volb,
		     const Matrix& indices);
  template<class T>
  void indexmultiply(volume<T>& vola, const volume<T>& volb,
		     const Matrix& indices);
  template<class T>
  void indexdivide(volume<T>& vola, const volume<T>& volb,
		   const Matrix& indices);
  
  template<class T>
  void indexset(volume<T>& vola, const Matrix& indices,
		const T num);

  template<class T>
  void indexadd(volume4D<T>& vola, const volume4D<T>& volb,
		const Matrix& indices);
  template<class T>
  void indexsubtract(volume4D<T>& vola, const volume4D<T>& volb,
		     const Matrix& indices);
  template<class T>
  void indexmultiply(volume4D<T>& vola, const volume4D<T>& volb,
		     const Matrix& indices);
  template<class T>
  void indexdivide(volume4D<T>& vola, const volume4D<T>& volb,
		   const Matrix& indices);
  
  template<class T>
  void indexset(volume4D<T>& vola, const Matrix& indices, 
	        const T num);

  
  // General Mathematical Operations

  template <class T>
  void clamp(volume<T>& vol, T minval, T maxval);
  template <class T>
  void clamp(volume4D<T>& vol, T minval, T maxval);

  template <class T>
  volume<T> binarise(const volume<T>& vol, T lowerth, T upperth, threshtype tt=inclusive);
  template <class T>
  volume<T> binarise(const volume<T>& vol, T thresh);
  template <class T>
  volume4D<T> binarise(const volume4D<T>& vol, T lowerth, T upperth, threshtype tt=inclusive);
  template <class T>
  volume4D<T> binarise(const volume4D<T>& vol, T thresh);

  template <class T>
  volume<T> threshold(const volume<T>& vol, T lowerth, T upperth, threshtype tt=inclusive);
  template <class T>
  volume<T> threshold(const volume<T>& vol, T thresh);
  template <class T>
  volume4D<T> threshold(const volume4D<T>& vol, T lowerth, T upperth, threshtype tt=inclusive);
  template <class T>
  volume4D<T> threshold(const volume4D<T>& vol, T thresh);

  template <class T>
  volume<T> min(const volume<T>& v1, const volume<T>& v2);
  template <class T>
  volume<T> max(const volume<T>& v1, const volume<T>& v2);
  template <class T>
  volume4D<T> min(const volume4D<T>& v1, const volume4D<T>& v2);
  template <class T>
  volume4D<T> max(const volume4D<T>& v1, const volume4D<T>& v2);

  volume<float> sqrt(const volume<char>& vol);
  volume<float> sqrt(const volume<short>& vol);
  volume<float> sqrt(const volume<int>& vol);
  volume<float> sqrt(const volume<float>& vol);
  volume<double> sqrt(const volume<double>& vol);
  template <class T>
  volume<float> sqrt_float(const volume<T>& vol);
 
  volume4D<float> sqrt(const volume4D<char>& vol4);
  volume4D<float> sqrt(const volume4D<short>& vol4);
  volume4D<float> sqrt(const volume4D<int>& vol4);
  volume4D<float> sqrt(const volume4D<float>& vol4);
  volume4D<double> sqrt(const volume4D<double>& vol4);
  template <class T>
  volume4D<float> sqrt_float(const volume4D<T>& vol4);

  template <class T>
  volume<float> meanvol(const volume4D<T>& vol4);
  template <class T>
  volume<float> stddevvol(const volume4D<T>& vol4);
  template <class T>
  volume<float> variancevol(const volume4D<T>& vol4);
  template <class T>
  volume<float> sumvol(const volume4D<T>& vol4);
  template <class T>
  volume<float> sumsquaresvol(const volume4D<T>& vol4);
  template <class T>
  volume<float> dotproductvol(const volume4D<T>& vol4, 
			      const ColumnVector& vec);

  template <class T>
  void pad(const volume<T>& vol, volume<T>& paddedvol);
  template <class T>
  void pad(const volume<T>& vol, volume<T>& paddedvol, 
	     int offsetx, int offsety, int offsetz);

  // Considers each volume as a vector and returns the dotproduct
  template <class T>
  double dotproduct(const volume<T>&  vol1,
                    const volume<T>&  vol2);
  template <class T, class S>
  double dotproduct(const volume<T>&  vol1,
                    const volume<T>&  vol2,
                    const volume<S>&  mask);
  template <class T, class S>
  double dotproduct(const volume<T>& vol1,
                    const volume<T>& vol2,
                    const volume<S>  *mask);

  // This is a bit of a special-needs funtion. If we consider vol1 and vol2
  // as column vectors v1 and v2 then powerdotproduct returns (v1.^n1)' * (v2.^n2)
  template <class T>
  double powerdotproduct(const volume<T>&  vol1,
                         unsigned int      n1,
                         const volume<T>&  vol2,
                         unsigned int      n2);  
  template <class T, class S>
  double powerdotproduct(const volume<T>&  vol1,
                         unsigned int      n1,
                         const volume<T>&  vol2,
                         unsigned int      n2,
                         const volume<S>&  mask);  
  template <class T, class S>
  double powerdotproduct(const volume<T>&  vol1,
                         unsigned int      n1,
                         const volume<T>&  vol2,
                         unsigned int      n2,
                         const volume<S>   *mask);  

  // These are global functions that duplicate some of the member functions.
  // The reason is that we want to be able to double template them in a 
  // convenient manner (e.g. use a <char> mask for <float> data).
  template <class T>
  double mean(const volume<T>&  vol);
  template <class T, class S>
  double mean(const volume<T>&  vol,
              const volume<S>&  mask);
  template <class T, class S>
  double mean(const volume<T>& vol,
              const volume<S> *mask);

  template <class T>
  double sum(const volume<T>&  vol);
  template <class T, class S>
  double sum(const volume<T>&  vol,
             const volume<S>&  mask);
  template <class T, class S>
  double sum(const volume<T>& vol,
             const volume<S> *mask);


  ///////////////////////////////////////////////////////////////////////////

  // IMAGE PROCESSING ROUTINES

  template <class T>
  void affine_transform(const volume<T>& vin, volume<T>& vout,
			     const Matrix& aff, float paddingsize=0.0);
  template <class T>
  void affine_transform(const volume<T>& vin, volume<T>& vout,
			const Matrix& aff, interpolation interptype, 
			float paddingsize=0.0);
  template <class T>
  volume<T> affine_transform_mask(const volume<T>& vin, const volume<T>& vout,
				  const Matrix& aff, float padding=0.0);

  template <class T>
  void get_axis_orientations(const volume<T>& inp1, 
			     int& icode, int& jcode, int& kcode);


  // the following convolve function do not attempt to normalise the kernel
  template <class T, class S>
  volume<T> convolve(const volume<T>& source, const volume<S>& kernel);
  template <class T, class S>
  volume<T> convolve_separable(const volume<T>& source, 
			       const ColumnVector& kernelx, 
			       const ColumnVector& kernely,
			       const ColumnVector& kernelz);

  // the following convolve functions take in a mask and also renormalise
  //  the result according to the overlap of kernel and mask at each point
  // NB: these functions should NOT be used with zero-sum kernels (eg.Laplacian)
  template <class T, class S, class M>
  volume<T> convolve(const volume<T>& source, const volume<S>& kernel, 
		     const volume<M>& mask, bool ignoremask=false, bool renormalise=true);
  template <class T, class M>
  volume<T> convolve_separable(const volume<T>& source, 
			       const ColumnVector& kernelx, 
			       const ColumnVector& kernely,
			       const ColumnVector& kernelz,
			       const volume<M>& mask, bool ignoremask=false, bool renormalise=true);

  //This implements Professor Smith's SUSAN convolve algorithm, note the number of optional parameters
  template <class T, class S>
    volume<T> susan_convolve(volume<T> source, const volume<S>& kernel, const float sigmabsq, const bool use_median, int num_usan,volume<T>* usan_area = new volume<T>(1,1,1),volume<T> usan_vol1=volume<T>(1,1,1),const float sigmab1sq=0,volume<T> usan_vol2 = volume<T>(1,1,1),const float sigmab2sq=0);

  template <class T, class S> 
    volume4D<T> generic_convolve(const volume4D<T>& source, const volume<S>& kernel, bool seperable=false, bool renormalise=true);
  template <class T, class S>
    volume<T> efficient_convolve(const volume<T>& vin, const volume<S>& vker);
  template <class T, class S>
  int insertpart(volume<T>& v1, const volume<S>& v2); 
  template <class T, class S, class U>
  volume<S> extractpart(const volume<T>& v1, const volume<S>& v2, const volume<U>& kernel) ;
   float fsllog2(float x);

    
   template <class T>
     volume4D<T> bandpass_temporal_filter(volume4D<T>& source,double hp_sigma, double lp_sigma);


  template <class T, class S>
  volume<T> morphfilter(const volume<T>& source, const volume<S>& mask,
			const string& filtertype);

  template <class T>
  volume<T> isotropic_resample(const volume<T>& aniso, float scale);

  template <class T>
  volume<T> subsample_by_2(const volume<T>& refvol, bool centred=true);
  template <class T>
  volume4D<T> subsample_by_2(const volume4D<T>& refvol, bool centred=true);

  void make_blur_mask(ColumnVector& bmask, const float final_vox_dim, 
		     const float init_vox_dim);
  template <class T>
  volume<T> blur(const volume<T>& source, const ColumnVector& resel_size);
  template <class T>
  volume<T> blur(const volume<T>& source, float iso_resel_size);
  template <class T>
  volume4D<T> blur(const volume4D<T>& source, const ColumnVector& resel_size);
  template <class T>
  volume4D<T> blur(const volume4D<T>& source, float iso_resel_size);
  template <class T>
  volume<T> smooth(const volume<T>& source, float sigma_mm);
  template <class T>
  volume<T> smooth2D(const volume<T>& source, float sigma_mm, int nulldir=3);
  template <class T>
  volume4D<T> smooth(const volume4D<T>& source, float sigma_mm);
  template <class T>
  volume4D<T> smooth2D(const volume4D<T>& source, float sigma_mm, int nulldir=3);

  ColumnVector gaussian_kernel1D(float sigma, int radius);
  volume<float> gaussian_kernel2D(float sigma, int radius);
  volume<float> gaussian_kernel3D(float sigma, int radius);
  volume<float> gaussian_kernel3D(float sigma,float xdim,float ydim,float zdim,float cutoff=4.0);
  volume<float> spherical_kernel(float radius, float xdim, float ydim, float zdim);
  volume<float> box_kernel(float length,float xdim,float ydim,float zdim);  //mm dimensions
  volume<float> box_kernel(int x,int y, int z);                        //voxel dimensions

  void make_grad_masks(volume<float>& maskx, volume<float>& masky, 
		       volume<float>& maskz);

  template <class T>
  volume<float> gradient(const volume<T>& source);

  template <class T>
  void gradient(const volume<T>& source,volume4D<float>& grad);

  // separate left and right gradients (changes at voxel mid-point)
  template <class T>
  volume4D<float> lrxgrad(const volume<float>& im, const volume<T>& mask);
  template <class T>
  volume4D<float> lrygrad(const volume<float>& im, const volume<T>& mask);
  template <class T>
  volume4D<float> lrzgrad(const volume<float>& im, const volume<T>& mask);
  
  
  template <class T>
  volume<T> log_edge_detect(const volume<T>& source, 
			    float sigma1, float sigma2, 
			    float zero_tolerance, bool twodimensional=false);
  template <class T>
  volume<T> fixed_edge_detect(const volume<T>& source, float threshold, 
			      bool twodimensional=false);
  template <class T>
  volume4D<T> edge_strengthen(const volume4D<T>& source);



  template <class T>
    volume<int> connected_components(const volume<T>& vol,  ColumnVector& clustersize, int numconnected=26);
  
  template <class T>
    volume<int> connected_components(const volume<T>& vol, 
                                   const volume<T>& mask, 
                                   bool (*binaryrelation)(T , T), ColumnVector& clustersize);

  template <class T>
  volume<int> connected_components(const volume<T>& vol, 
				   int numconnected=26);
  template <class T>
  volume<int> connected_components(const volume<T>& vol, 
                                   const volume<T>& mask, 
                                   bool (*binaryrelation)(T , T));

  template <class T>
  volume<float> distancemap(const volume<T>& binaryvol);
  template <class T>
  volume<float> distancemap(const volume<T>& binaryvol, const volume<T>& mask);

  template <class T>
  volume4D<float> sparseinterpolate(const volume4D<T>& sparsesamps, const volume<T>& mask,
				    const string& interpmethod="general");
  // can have "general" or "nearestneighbour" (or "nn") for interpmethod


 template <class T>
 Matrix NewimageVox2NewimageVoxMatrix(const Matrix& flirt_in2ref,
		      const volume<T>& invol, const volume<T>& refvol);



  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////



  ///////////////////////////////////////////////////////////////////////////
  // TEMPLATE DEFINITIONS
  ///////////////////////////////////////////////////////////////////////////

  template <class T>
  void print_size(const volume<T>& source)
    {
      cout << "Size: " << source.xsize() << ", " << source.ysize()
	   << ", " << source.zsize() << endl;
    }


  template <class T>
  void print_volume_info(const volume<T>& source, const string& name)
    {
      cout << name << ":: Size = (" << source.xsize() << ","
	   << source.ysize() << "," << source.zsize() << ")" << endl;
      cout << name << ":: ROI Size = (" << source.maxx() - source.minx() + 1 
	   << "," << source.maxy() - source.miny() + 1
	   << "," << source.maxz() - source.minz() + 1 << ")" << endl;
      cout << name << ":: Dims = (" << source.xdim() << ","
	   << source.ydim() << "," << source.zdim() << ")" << endl;
      cout << name << ":: Minimum and maximum intensities are: " 
	   << source.min() << " and " << source.max() << endl;
    }  

  template <class T>
  void print_volume_info(const volume4D<T>& source, const string& name)
    {
      cout << name << ":: Size = (" << source.xsize() << ","
	   << source.ysize() << "," << source.zsize() << ","
	   << source.tsize() << ")" << endl;
      cout << name << ":: ROI Size = (" << source.maxx() - source.minx() + 1 
	   << "," << source.maxy() - source.miny() + 1
	   << "," << source.maxz() - source.minz() + 1 
	   << "," << source.maxt() - source.mint() + 1 << ")" << endl;
      cout << name << ":: Dims = (" << source.xdim() << ","
	   << source.ydim() << "," << source.zdim() << "," 
	   << source.tdim() << ")" << endl;
      cout << name << ":: Minimum and maximum intensities are: " 
	   << source.min() << " and " << source.max() << endl;
    }  

  
  ///////////////////////////////////////////////////////////////////////////
  // BASIC IMAGE SUPPORT FUNCTIONS
  ///////////////////////////////////////////////////////////////////////////

  template <class T>
  void clamp(volume<T>& vol, T minval, T maxval)
    {
      if (maxval < minval) return;
      for (int z=vol.minz(); z<=vol.maxz(); z++) {
	for (int y=vol.miny(); y<=vol.maxy(); y++) {
	  for (int x=vol.minx(); x<=vol.maxx(); x++) {
	    if (vol.value(x,y,z)>maxval)       vol.value(x,y,z)=maxval;
	    else if (vol.value(x,y,z)<minval)  vol.value(x,y,z)=minval;
	  }
	}
      }
    }
  
  template <class T>
  void clamp(volume4D<T>& vol, T minval, T maxval) 
  {
    for (int t=vol.mint(); t<=vol.maxt(); t++) {
      clamp(vol[t],minval,maxval);
    }
  }


  template <class T>
  volume<T> abs(const volume<T>& vol)
    {
      volume<T> newvol(vol);
      for (int z=newvol.minz(); z<=newvol.maxz(); z++) {
	for (int y=newvol.miny(); y<=newvol.maxy(); y++) {
	  for (int x=newvol.minx(); x<=newvol.maxx(); x++) {
	    newvol.value(x,y,z) = (T) fabs((double) newvol.value(x,y,z));
	  }
	}
      }
      return newvol;
    }

  template <class T>
  volume4D<T> abs(const volume4D<T>& vol)
    {
      volume4D<T> newvol(vol);
      for (int z=newvol.minz(); z<=newvol.maxz(); z++) {
	for (int y=newvol.miny(); y<=newvol.maxy(); y++) {
	  for (int x=newvol.minx(); x<=newvol.maxx(); x++) {
	    for (int t=newvol.mint(); t<=newvol.maxt(); t++) {
	      newvol(x,y,z,t) = (T) fabs((double) newvol.value(x,y,z,t));
	    }
	  }
	}
      }
      return newvol;
    }

  template <class T>
  volume<T> binarise(const volume<T>& vol, T lowerth, T upperth, threshtype tt)
    { 
      volume<T> newvol(vol);
      newvol.binarise(lowerth,upperth,tt);
      return newvol;
    }

  template <class T>
  volume<T> binarise(const volume<T>& vol, T lowthresh)
    { 
      return binarise(vol,lowthresh,vol.max(),inclusive); 
    }

  template <class T>
  volume4D<T> binarise(const volume4D<T>& vol, T lowerth, T upperth, threshtype tt)
    { 
      volume4D<T> newvol(vol);
      newvol.binarise(lowerth,upperth,tt);
      return newvol;
    }

  template <class T>
  volume4D<T> binarise(const volume4D<T>& vol, T thresh)
    { 
      return binarise(vol,thresh,vol.max(),inclusive); 
    }


  template <class T>
  volume<T> threshold(const volume<T>& vol, T lowerth, T upperth, threshtype tt)
    { 
      volume<T> newvol(vol);
      newvol.threshold(lowerth,upperth,tt);
      return newvol;
    }

  template <class T>
  volume<T> threshold(const volume<T>& vol, T thresh)
    { 
      return threshold(vol,thresh,vol.max(),inclusive); 
    }

  template <class T>
  volume4D<T> threshold(const volume4D<T>& vol, T lowerth, T upperth, threshtype tt)
    { 
      volume4D<T> newvol(vol);
      newvol.threshold(lowerth,upperth,tt);
      return newvol;
    }

  template <class T>
  volume4D<T> threshold(const volume4D<T>& vol, T thresh)
    { 
      return threshold(vol,thresh,vol.max(),inclusive); 
    }

  template <class T>
  volume<T> min(const volume<T>& v1, const volume<T>& v2) 
    {
      if (!samesize(v1,v2)) {
	imthrow("Must use volumes of same size in min(v1,v2)",3);
      }
      volume<T> newvol(v1);
      for (int z=newvol.minz(); z<=newvol.maxz(); z++) {
	for (int y=newvol.miny(); y<=newvol.maxy(); y++) {
	  for (int x=newvol.minx(); x<=newvol.maxx(); x++) {
	    newvol(x,y,z) = Min(v1(x,y,z),v2(x,y,z));
	  }
	}
      }
      return newvol;
    }

  template <class T>
  volume<T> max(const volume<T>& v1, const volume<T>& v2)
    {
      if (!samesize(v1,v2)) {
	imthrow("Must use volumes of same size in min(v1,v2)",3);
      }
      volume<T> newvol(v1);
      for (int z=newvol.minz(); z<=newvol.maxz(); z++) {
	for (int y=newvol.miny(); y<=newvol.maxy(); y++) {
	  for (int x=newvol.minx(); x<=newvol.maxx(); x++) {
	    newvol(x,y,z) = Max(v1(x,y,z),v2(x,y,z));
	  }
	}
      }
      return newvol;
    }

  template <class T>
  volume4D<T> min(const volume4D<T>& v1, const volume4D<T>& v2)
    {
      if (!samesize(v1,v2)) {
	imthrow("Must use volumes of same size in min(v1,v2)",3);
      }
      volume4D<T> newvol(v1);
      for (int t=newvol.mint(); t<=newvol.maxt(); t++) {
	newvol[t] = min(v1[t],v2[t]);
      }
      return newvol;
    }

  template <class T>
  volume4D<T> max(const volume4D<T>& v1, const volume4D<T>& v2)
    {
      if (!samesize(v1,v2)) {
	imthrow("Must use volumes of same size in min(v1,v2)",3);
      }
      volume4D<T> newvol(v1);
      for (int t=newvol.mint(); t<=newvol.maxt(); t++) {
	newvol[t] = max(v1[t],v2[t]);
      }
      return newvol;
    }


  template <class T>
  volume<float> sqrt_float(const volume<T>& vol)
  {
    volume<float> retvol;
    copyconvert(vol,retvol);
    for (int z=vol.minz(); z<=vol.maxz(); z++) {
      for (int y=vol.miny(); y<=vol.maxy(); y++) {
	for (int x=vol.minx(); x<=vol.maxx(); x++) {
	  if (vol(x,y,z)>0) {
	    retvol(x,y,z) = sqrt((double) vol(x,y,z));
	  } else {
	    retvol(x,y,z) = 0;
	  }
	}
      }
    }
    return retvol;
  }


  template <class T>
  volume4D<float> sqrt_float(const volume4D<T>& vol4)
  {
    if (vol4.mint()<0) { volume4D<float> newvol; return newvol; }
    volume4D<float> retvol;
    copyconvert(vol4,retvol);
    for (int t=vol4.mint(); t<=vol4.maxt(); t++) {
      for (int z=vol4.minz(); z<=vol4.maxz(); z++) {
	for (int y=vol4.miny(); y<=vol4.maxy(); y++) {
	  for (int x=vol4.minx(); x<=vol4.maxx(); x++) {
	    if (vol4(x,y,z,t)>0) {
	      retvol(x,y,z,t) = sqrt((double) vol4(x,y,z,t));
	    } else {
	      retvol(x,y,z,t) = 0;
	    }
	  }
	}
      }
    }
    return retvol;
  }


  template <class T>

  volume<float> sumvol(const volume4D<T>& vol4)
  {
    if (vol4.mint()<0) { volume<float> newvol; return newvol; }
    volume<float> SumVol, dummy;
    copyconvert(vol4[vol4.mint()],SumVol);
    for (int ctr=vol4.mint() + 1; ctr <= vol4.maxt(); ctr++) {
      copyconvert(vol4[ctr],dummy);
      SumVol += dummy;
    }       
    return SumVol;
  }


  template <class T>
  volume<float> meanvol(const volume4D<T>& vol4)
  {
    if (vol4.mint()<0) { volume<float> newvol; return newvol; }
    volume<float> MeanVol = sumvol(vol4);
    if (vol4.maxt() > vol4.mint()) 
      MeanVol /= (float) (vol4.maxt() - vol4.mint()+ 1);
    return MeanVol;
  }


  template <class T>
  volume<float> sumsquaresvol(const volume4D<T>& vol4)
  {
    if (vol4.mint()<0) { volume<float> newvol; return newvol; }
    volume<float> SumSq, dummy;
    copyconvert(vol4[vol4.mint()] * vol4[vol4.mint()],SumSq);
    for (int ctr=vol4.mint() + 1; ctr <= vol4.maxt(); ctr++) {
      copyconvert(vol4[ctr] * vol4[ctr],dummy);
      SumSq += dummy;
    }       
    return SumSq;
  }
  //This rewrite of variancevol gives the same output as doing the sumsquaresvol etc in double precision internally
  template <class T>
  volume<float> variancevol(const volume4D<T>& vol4)
  {
     volume<float> variance;
     if (vol4.mint()<0) 
       return variance; 
     volume<float> Mean = meanvol(vol4);
     variance.reinitialize(vol4.xsize(),vol4.ysize(),vol4.zsize());
     float n=1.0;
     if (vol4.maxt() > vol4.mint()) { n = (float) (vol4.maxt() - vol4.mint() + 1); }

     for (int z=vol4.minz(); z<=vol4.maxz(); z++) 
       for (int y=vol4.miny(); y<=vol4.maxy(); y++) 
	 for (int x=vol4.minx(); x<=vol4.maxx(); x++) { 
	   double total(0);
	   for (int t=vol4.mint(); t<=vol4.maxt(); t++) 
	     total+=pow(vol4(x,y,z,t)-Mean(x,y,z),2.0);
	   variance(x,y,z)=(float)total;
	 }
     variance /= (float) (n-1.0);
     return variance;
  }

  template <class T>
  volume<float> stddevvol(const volume4D<T>& vol4)
  {
    if (vol4.mint()<0) { volume<float> newvol; return newvol; }
    volume<float> StdVol = variancevol(vol4);
    for (int z=StdVol.minz(); z<=StdVol.maxz(); z++) {
      for (int y=StdVol.miny(); y<=StdVol.maxy(); y++) {
	for (int x=StdVol.minx(); x<=StdVol.maxx(); x++) {
	  StdVol(x,y,z) = sqrt(StdVol(x,y,z));
	}
      }
    }
    return StdVol;
  }

  template <class T>
  volume<float> dotproductvol(const volume4D<T>& vol4, 
			      const ColumnVector& vec)
  {
    if (vol4.mint()<0) { volume<float> newvol; return newvol; }
    if ( (vol4.maxt() - vol4.mint() + 1) != vec.Nrows() )
      {
	cerr << "ERROR::Time series length differs from vector length in"
	     << " dotproductvol()" << endl;
	volume<float> newvol; return newvol;
      }
    volume<float> vol4copy;
    copyconvert(vol4[vol4.mint()],vol4copy);
    volume<float> DotVol(vol4copy);
    DotVol *= (float) vec(1);
    for (int n=vol4.mint() + 1; n <= vol4.maxt(); n++) {
      copyconvert(vol4[n],vol4copy);
      DotVol += (vol4copy * (float) vec(1 + n - vol4.mint()));
    }
    return DotVol;
  }


  template <class T>
  void pad(const volume<T>& vol, volume<T>& paddedvol)
    {
      // The default type of padding is central padding
      int wx1, wy1, wz1, wx2, wy2, wz2, offx, offy, offz;
      wx1 = vol.maxx() - vol.minx() + 1;
      wy1 = vol.maxy() - vol.miny() + 1;
      wz1 = vol.maxz() - vol.minz() + 1;
      wx2 = paddedvol.maxx() - paddedvol.minx() + 1;
      wy2 = paddedvol.maxy() - paddedvol.miny() + 1;
      wz2 = paddedvol.maxz() - paddedvol.minz() + 1;
      if ( (wx2<wx1) || (wy2<wy1) || (wz2<wz1) ) {
	imthrow("Cannot pad when target volume is smaller than original",7);
      }
      offx = (wx2 - wx1)/2;
      offy = (wy2 - wy1)/2;
      offz = (wz2 - wz1)/2;
      pad(vol,paddedvol,offx,offy,offz);
    }

  template <class T>
  double dotproduct(const volume<T>& vol1,
                    const volume<T>& vol2)
  {
    if (!samesize(vol1,vol2,true)) imthrow("dotproduct: Image dimension mismatch",99);

    double rval = 0.0;
    for (typename volume<T>::fast_const_iterator it1=vol1.fbegin(), it_end=vol1.fend(), it2=vol2.fbegin(); it1 != it_end; ++it1, ++it2) {
      rval += static_cast<double>((*it1)*(*it2));
    }
    return(rval);
  }
  
  template <class T, class S>
  double dotproduct(const volume<T>& vol1,
                    const volume<T>& vol2,
                    const volume<S>& mask)
  {
    if (!samesize(vol1,vol2,true) || !samesize(vol1,mask,true)) imthrow("dotproduct: Image dimension mismatch",99);

    double rval = 0.0;
    typename volume<S>::fast_const_iterator  itm = mask.fbegin();
    for (typename volume<T>::fast_const_iterator it1=vol1.fbegin(), it_end=vol1.fend(), it2=vol2.fbegin(); it1 != it_end; ++it1, ++it2, ++itm) {
      if (*itm > 0.5) {
        rval += static_cast<double>((*it1)*(*it2));
      }
    }
    return(rval);
  }

  template <class T, class S>
  double dotproduct(const volume<T>& vol1,
                    const volume<T>& vol2,
                    const volume<S>  *mask)
  {
    if (mask) return(dotproduct(vol1,vol2,*mask));
    else return(dotproduct(vol1,vol2));
  }

  template <class T>
  double powerdotproduct(const volume<T>&  vol1,
                         unsigned int      n1,
                         const volume<T>&  vol2,
                         unsigned int      n2)
  {
    if (!samesize(vol1,vol2,true)) imthrow("powerdotproduct: Image dimension mismatch",99);

    double rval = 0.0;
    for (typename volume<T>::fast_const_iterator it1=vol1.fbegin(), it_end=vol1.fend(), it2=vol2.fbegin(); it1 != it_end; ++it1, ++it2) {
      double val1 = 1.0;
      for (unsigned int i=0; i<n1; i++) val1 *= *it1;
      double val2 = 1.0;
      for (unsigned int i=0; i<n2; i++) val2 *= *it2;
      rval += static_cast<double>(val1*val2);
    }
    return(rval);
  } 

  template <class T, class S>
  double powerdotproduct(const volume<T>&  vol1,
                         unsigned int      n1,
                         const volume<T>&  vol2,
                         unsigned int      n2,
                         const volume<S>&  mask)
  {
    if (!samesize(vol1,vol2,true) || !samesize(vol1,mask,true)) imthrow("powerdotproduct: Image dimension mismatch",99);

    double rval = 0.0;
    typename volume<S>::fast_const_iterator itm = mask.fbegin();
    for (typename volume<T>::fast_const_iterator it1=vol1.fbegin(), it_end=vol1.fend(), it2=vol2.fbegin(); it1 != it_end; ++it1, ++it2, ++itm) {
      if (*itm > 0.5) {
        double val1 = 1.0;
        for (unsigned int i=0; i<n1; i++) val1 *= *it1;
        double val2 = 1.0;
        for (unsigned int i=0; i<n2; i++) val2 *= *it2;
        rval += static_cast<double>(val1*val2);
      }
    }
    return(rval);
  } 

  template <class T, class S>
  double powerdotproduct(const volume<T>&  vol1,
                         unsigned int      n1,
                         const volume<T>&  vol2,
                         unsigned int      n2,
                         const volume<S>   *mask)
  {
    if (mask) return(powerdotproduct(vol1,n1,vol2,n2,*mask));
    else return(powerdotproduct(vol1,n1,vol2,n2));
  }

  template <class T>
  double mean(const volume<T>& vol)
  {
    double rval = sum(vol);
    rval /= static_cast<double>(vol.xsize()*vol.ysize()*vol.zsize());
    return(rval);
  }

  template <class T, class S>
  double mean(const volume<T>& vol,
              const volume<S>& mask)
  {
    if (!samesize(vol,mask,true)) imthrow("mean: Image-Mask dimension mismatch",99);

    double rval = 0.0;
    unsigned int n = 0;
    typename volume<S>::fast_const_iterator  itm = mask.fbegin();
    for (typename volume<T>::fast_const_iterator it=vol.fbegin(), it_end=vol.fend(); it != it_end; ++it, ++itm) {
      if (*itm > 0.5) {
        n++;
        rval += static_cast<double>(*it);
      }
    }
    rval /= static_cast<double>(n);
    return(rval);
  }

  template <class T, class S>
  double mean(const volume<T>& vol,
              const volume<S> *mask)
  {
    if (mask) return(mean(vol,*mask));
    else return(mean(vol));
  }

  template <class T>
  double sum(const volume<T>& vol)
  {
    double rval = 0.0;
    for (typename volume<T>::fast_const_iterator it=vol.fbegin(), it_end=vol.fend(); it != it_end; ++it) {
      rval += static_cast<double>(*it);
    }
    return(rval);
  }

  template <class T, class S>
  double sum(const volume<T>& vol,
             const volume<S>& mask)
  {
    if (!samesize(vol,mask,true)) imthrow("sum: Image-Mask dimension mismatch",99);

    double rval = 0.0;
    typename volume<S>::fast_const_iterator  itm = mask.fbegin();
    for (typename volume<T>::fast_const_iterator it=vol.fbegin(), it_end=vol.fend(); it != it_end; ++it, ++itm) {
      if (*itm > 0.5) {
        rval += static_cast<double>(*it);
      }
    }
    return(rval);
  }

  template <class T, class S>
  double sum(const volume<T>& vol,
              const volume<S> *mask)
  {
    if (mask) return(sum(vol,*mask));
    else return(sum(vol));
  }

  template <class T, class S>
  volume<T> divide(const volume<T>& numervol, const volume<T>& denomvol,
		   const volume<S>& mask)
  {
    if ((!samesize(numervol,denomvol)) || (!samesize(mask,denomvol))) {
      imthrow("Attempted to divide images/ROIs of different sizes",3);
    }
    volume<T> resvol(numervol);
    for (int z=denomvol.minz(); z<=denomvol.maxz(); z++) {
      for (int y=denomvol.miny(); y<=denomvol.maxy(); y++) {
	for (int x=denomvol.minx(); x<=denomvol.maxx(); x++) {
	  if (mask(x,y,z)!=0) {
	    resvol(x,y,z) /= denomvol(x,y,z);
	  } else { 
	    resvol(x,y,z) = 0;
	  }
	}
      }
    }
    return resvol;
  }


  template <class T, class S>
  volume4D<T> divide(const volume4D<T>& numervol, const volume4D<T>& denomvol,
		   const volume<S>& mask)
  {
    if ( (denomvol.tsize()<1) ||
	 (!samesize(numervol,denomvol)) || (!samesize(mask,denomvol[0]))) {
      imthrow("Attempted to divide images/ROIs of different sizes",3);
    }
    volume4D<T> resvol(numervol);
    for (int z=denomvol.minz(); z<=denomvol.maxz(); z++) {
      for (int y=denomvol.miny(); y<=denomvol.maxy(); y++) {
	for (int x=denomvol.minx(); x<=denomvol.maxx(); x++) {
	  if (mask(x,y,z)!=0) {
	    for (int t=denomvol.mint(); t<=denomvol.maxt(); t++) {
	      resvol(x,y,z,t) /= denomvol(x,y,z,t);
	    }
	  } else {
	    for (int t=denomvol.mint(); t<=denomvol.maxt(); t++) {
	      resvol(x,y,z,t) = 0;
	    }
	  }
	}
      }
    }
    return resvol;
  }


  template <class T, class S>
  volume4D<T> divide(const volume4D<T>& numervol, const volume<T>& denomvol,
		   const volume<S>& mask)
  {
    if ( (numervol.tsize()<1) || 
         (!samesize(numervol[0],denomvol)) || (!samesize(mask,denomvol))) {
      imthrow("Attempted to divide images/ROIs of different sizes",3);
    }
    volume4D<T> resvol(numervol);
    for (int z=denomvol.minz(); z<=denomvol.maxz(); z++) {
      for (int y=denomvol.miny(); y<=denomvol.maxy(); y++) {
	for (int x=denomvol.minx(); x<=denomvol.maxx(); x++) {
	  if (mask(x,y,z)!=0) {
	    for (int t=numervol.mint(); t<=numervol.maxt(); t++) {
	      resvol(x,y,z,t) /= denomvol(x,y,z);
	    }
	  } else {
	    for (int t=numervol.mint(); t<=numervol.maxt(); t++) {
	      resvol(x,y,z,t) = 0;
	    }
	  }	    
	}
      }
    }
    return resvol;
  }

template <class T, class S>
  volume<T> mask_volume( const volume<T>& invol, const volume<S>& mask)
  {
    if ( !samesize(invol,mask)) {
      imthrow("Attempted to mask with wrong sized mask",3);
    }
    volume<T> resvol(invol);
    for (int z=invol.minz(); z<=invol.maxz(); z++) {
      for (int y=invol.miny(); y<=invol.maxy(); y++) {
	for (int x=invol.minx(); x<=invol.maxx(); x++) {
	  if (mask(x,y,z)!=0.0) {
	    resvol(x,y,z)=invol(x,y,z);
	  } else {
	    resvol(x,y,z) = (T)0.0;   
	  }	    
	}
      }
    }
    return resvol;
  }



  template<class T>
  void indexadd(volume<T>& vola, const volume<T>& volb, const Matrix& indices)
  {
    if ((!samesize(vola,volb))) {
      imthrow("Attempted to add images/ROIs of different sizes",3);
    }
    if (indices.Ncols()!=3) {
      imthrow("indexadd: must specify indices in an Nx3 matrix",11);
    }
    for (int n=1; n<=indices.Nrows(); n++) {
      vola((int)indices(n,1),(int)indices(n,2),(int)indices(n,3)) += 
	volb((int)indices(n,1),(int)indices(n,2),(int)indices(n,3));
    }
  }


  template<class T>
  void indexsubtract(volume<T>& vola, const volume<T>& volb,
		     const Matrix& indices)
  {
    if ((!samesize(vola,volb))) {
      imthrow("Attempted to subtract images/ROIs of different sizes",3);
    }
    if (indices.Ncols()!=3) {
      imthrow("indexsubtract: must specify indices in an Nx3 matrix",11);
    }
    for (int n=1; n<=indices.Nrows(); n++) {
      vola((int)indices(n,1),(int)indices(n,2),(int)indices(n,3)) -= 
	volb((int)indices(n,1),(int)indices(n,2),(int)indices(n,3));
    }
  }

  template<class T>
  void indexmultiply(volume<T>& vola, const volume<T>& volb,
		     const Matrix& indices)
  {
    if ((!samesize(vola,volb))) {
      imthrow("Attempted to multiply images/ROIs of different sizes",3);
    }
    if (indices.Ncols()!=3) {
      imthrow("indexmultiply: must specify indices in an Nx3 matrix",11);
    }
    for (int n=1; n<=indices.Nrows(); n++) {
      vola((int)indices(n,1),(int)indices(n,2),(int)indices(n,3)) *= 
	volb((int)indices(n,1),(int)indices(n,2),(int)indices(n,3));
    }
  }

  template<class T>
  void indexdivide(volume<T>& vola, const volume<T>& volb,
		   const Matrix& indices)
  {
    if ((!samesize(vola,volb))) {
      imthrow("Attempted to divide images/ROIs of different sizes",3);
    }
    if (indices.Ncols()!=3) {
      imthrow("indexdivide: must specify indices in an Nx3 matrix",11);
    }
    for (int n=1; n<=indices.Nrows(); n++) {
      vola((int)indices(n,1),(int)indices(n,2),(int)indices(n,3)) /= 
	volb((int)indices(n,1),(int)indices(n,2),(int)indices(n,3));
    }
  }

  template<class T>
  void indexset(volume<T>& vola, const Matrix& indices, 
	        const T num)
  {

    if (indices.Ncols()!=3) {
      imthrow("indexset: must specify indices in an Nx3 matrix",11);
    }
    for (int n=1; n<=indices.Nrows(); n++) {
      vola((int)indices(n,1),(int)indices(n,2),(int)indices(n,3)) = num; 

    }
  }

  template<class T>
  void indexadd(volume4D<T>& vola, const volume4D<T>& volb,
		const Matrix& indices)
  {
    if ((!samesize(vola,volb))) {
      imthrow("Attempted to add images/ROIs of different sizes",3);
    }
    if ( (indices.Ncols()!=3) && (indices.Ncols()!=4) ) {
      imthrow("indexadd: must specify indices in an Nx3 or Nx4 matrix",11);
    }
    int tmax = vola.tsize(), tpt=0;
    if (indices.Ncols()==3) { tmax = 1; }
    for (int n=1; n<=indices.Nrows(); n++) {
      for (int t=0; t<tmax; t++) {
	if (indices.Ncols()==4) { tpt = (int)indices(n,4); } else { tpt = t; }
	vola((int)indices(n,1),(int)indices(n,2),(int)indices(n,3),tpt) += 
	  volb((int)indices(n,1),(int)indices(n,2),(int)indices(n,3),tpt);
      }
    }
  }

  template<class T>
  void indexsubtract(volume4D<T>& vola, const volume4D<T>& volb,
		     const Matrix& indices)
  {
    if ((!samesize(vola,volb))) {
      imthrow("Attempted to subtract images/ROIs of different sizes",3);
    }
    if ( (indices.Ncols()!=3) && (indices.Ncols()!=4) ) {
      imthrow("indexsubtract: must specify indices in an Nx3 or Nx4 matrix",11);
    }
    int tmax = vola.tsize(), tpt=0;
    if (indices.Ncols()==3) { tmax = 1; }
    for (int n=1; n<=indices.Nrows(); n++) {
      for (int t=0; t<tmax; t++) {
	if (indices.Ncols()==4) { tpt = (int)indices(n,4); } else { tpt = t; }
	vola((int)indices(n,1),(int)indices(n,2),(int)indices(n,3),tpt) -= 
	  volb((int)indices(n,1),(int)indices(n,2),(int)indices(n,3),tpt);
      }
    }
  }

  template<class T>
  void indexmultiply(volume4D<T>& vola, const volume4D<T>& volb,
		     const Matrix& indices)
  {
    if ((!samesize(vola,volb))) {
      imthrow("Attempted to multiply images/ROIs of different sizes",3);
    }
    if ( (indices.Ncols()!=3) && (indices.Ncols()!=4) ) {
      imthrow("indexmultiply: must specify indices in an Nx3 or Nx4 matrix",11);
    }
    int tmax = vola.tsize(), tpt=0;
    if (indices.Ncols()==3) { tmax = 1; }
    for (int n=1; n<=indices.Nrows(); n++) {
      for (int t=0; t<tmax; t++) {
	if (indices.Ncols()==4) { tpt = (int)indices(n,4); } else { tpt = t; }
	vola((int)indices(n,1),(int)indices(n,2),(int)indices(n,3),tpt) *= 
	  volb((int)indices(n,1),(int)indices(n,2),(int)indices(n,3),tpt);
      }
    }
  }

  template<class T>
  void indexdivide(volume4D<T>& vola, const volume4D<T>& volb,
		   const Matrix& indices)
  {
    if ((!samesize(vola,volb))) {
      imthrow("Attempted to divide images/ROIs of different sizes",3);
    }
    if ( (indices.Ncols()!=3) && (indices.Ncols()!=4) ) {
      imthrow("indexdivide: must specify indices in an Nx3 or Nx4 matrix",11);
    }
    int tmax = vola.tsize(), tpt=0;
    if (indices.Ncols()==3) { tmax = 1; }
    for (int n=1; n<=indices.Nrows(); n++) {
      for (int t=0; t<tmax; t++) {
	if (indices.Ncols()==4) { tpt = (int)indices(n,4); } else { tpt = t; }
	vola((int)indices(n,1),(int)indices(n,2),(int)indices(n,3),tpt) /= 
	  volb((int)indices(n,1),(int)indices(n,2),(int)indices(n,3),tpt);
      }
    }
  }

  template<class T>
  void indexset(volume4D<T>& vola, const Matrix& indices, const T num)
  {
    if ( (indices.Ncols()!=3) && (indices.Ncols()!=4) ) {
      imthrow("indexset: must specify indices in an Nx3 or Nx4 matrix",11);
    }
    int tmax = vola.tsize(), tpt=0;
    if (indices.Ncols()==3) { tmax = 1; }
    for (int n=1; n<=indices.Nrows(); n++) {
      for (int t=0; t<tmax; t++) {
	if (indices.Ncols()==4) { tpt = (int)indices(n,4); } else { tpt = t; }
	vola((int)indices(n,1),(int)indices(n,2),(int)indices(n,3),tpt) =0; 
      }
    }
  }

  // AFFINE TRANSFORM
 template <class T>
 void raw_affine_transform(const volume<T>& vin, volume<T>& vout,
			   const Matrix& aff);

  template <class T>
  void affine_transform_mask(const volume<T>& vin, volume<T>& vout,
			     const Matrix& aff, float padding, const T padval);

  template <class T>
  volume<T> affine_transform_mask(const volume<T>& vin, const volume<T>& vout,
				  const Matrix& aff, float padding)
    {
      volume<T> affmask;
      affmask = vout;
      affmask = (T) 1;
      affine_transform_mask(vin,affmask,aff,padding,(T) 0);
    }


  template <class T>
  void affine_transform(const volume<T>& vin, volume<T>& vout,
			const Matrix& aff, float paddingsize)
    {
      T padval = vin.getpadvalue();
      extrapolation oldex = vin.getextrapolationmethod();

      vin.setpadvalue(vin.backgroundval());
      vin.setextrapolationmethod(extraslice);

      raw_affine_transform(vin,vout,aff);
      // now mask the output to eliminate streaks formed by the sinc interp...
      affine_transform_mask(vin,vout,aff,paddingsize,vin.backgroundval());

      vin.setpadvalue(padval);
      vin.setextrapolationmethod(oldex);
    }

  template <class T>
  void affine_transform(const volume<T>& vin, volume<T>& vout,
			const Matrix& aff, interpolation interptype, 
			float paddingsize) 
    {
      interpolation oldinterp;
      oldinterp = vin.getinterpolationmethod();
      vin.setinterpolationmethod(interptype);
      affine_transform(vin,vout,aff,paddingsize);
      vin.setinterpolationmethod(oldinterp);
    }


  
  ///////////////////////////////////////////////////////////////////////////
  // CONVOLVE
    template <class T, class S>
    volume<T> susan_convolve(const volume<T> source, const volume<S>& kernel, const float sigmabsq, const bool use_median, int num_usan,volume<T>* usan_area ,volume<T> usan_vol1,const float sigmab1sq,volume<T> usan_vol2 ,const float sigmab2sq)
    //template <class T, class S, class U, class V, class W>
    //volume<T> susan_convolve(const volume<T>& source, const volume<S>& kernel, const float sigmabsq, const bool use_median, int num_usan,volume<U>* usan_area = new volume<T>(1,1,1),const volume<V>& usan_vol1=volume<T>(1,1,1),const float sigmab1sq=0,const volume<W>& usan_vol2 = volume<T>(1,1,1),const float sigmab2sq=0)
    //Note that the commented out declaration won't work with the optional arguements (since U,V,W need to be defined in call...). Code is provided for possible
    //future improvments, as it is all usans are templated as input
{
//need to use a pointer for usan_area as creating a default parameter for a pass-by-reference gives 
//a "assignment to tempory memory" warning in gcc
//default values for lut1 etc
  if (    (( (kernel.maxz() - kernel.minz()) % 2)==1) || 
	  (( (kernel.maxy() - kernel.miny()) % 2)==1) ||
	  (( (kernel.maxx() - kernel.minx()) % 2)==1) ) 
	  cerr << "WARNING:: Off-centre convolution being performed as kernel has even dimensions" << endl; 
  if ((num_usan>=1 && !samesize(source,usan_vol1)) || (num_usan>=2 && !samesize(source,usan_vol2)) || (num_usan>=1 && !samesize(source,*usan_area))      ) 
  {
    cerr << "Warning: an external usan or output usan is not the same size as the source image, reverting to num_usans=0 mode" << endl;
    num_usan=0;
  }
  int midx, midy, midz,lz,uz,lx,ux,ly,uy,lutsize=16384;
  lz=source.minz();
  uz=source.maxz();
  lx=source.minx();
  ux=source.maxx();
  ly=source.miny();
  uy=source.maxy();
  volume<T> result(source);
  midz=(kernel.maxz() - kernel.minz())/2;
  midy=(kernel.maxy() - kernel.miny())/2;
  midx=(kernel.maxx() - kernel.minx())/2;
  //generate look up table
  float range1=1,range2=1,range = (source.max() - source.min())/(float)lutsize;
  float **lut=new float *[3];
  for(int i=0;i<=num_usan;i++) lut[i]=new float[2*lutsize+1];
  for(int i=0;i<=num_usan;i++) lut[i]+=lutsize;
  for (int i=0;i<=lutsize;i++) lut[0][-i]=lut[0][i]= exp(-pow(i*range,2.0)/sigmabsq);
  if (num_usan>=1) 
  {
     range1= (usan_vol1.max() - usan_vol1.min())/(float)lutsize;
     for (int i=0;i<=lutsize;i++)   lut[1][-i]=lut[1][i]= exp(-pow(i*range1,2.0)/sigmab1sq);
  }
  if (num_usan>=2) 
  {
     range2= (usan_vol2.max() - usan_vol2.min())/(float)lutsize; 
     for (int i=0;i<=lutsize;i++)   lut[2][-i]=lut[2][i]= exp(-pow(i*range2,2.0)/sigmab2sq);
  }
  ColumnVector mediankernel((kernel.zsize()>1)?26:8); // cube or square, minus central voxel
  int medoffst=((kernel.zsize()>1)?1:0);
  for (int z=lz; z<=uz; z++) 
    for (int y=ly; y<=uy; y++) 
      for (int x=lx; x<=ux; x++) 
      {
	 int xmin=x-midx,ymin=y-midy,zmin=z-midz;
	 int xmax=MIN(x+midx,ux),ymax=MIN(y+midy,uy),zmax=MIN(z+midz,uz);
         float num=0, denom=0,center_val1=0,center_val2=0,factor;
         float center_val=source.value(x,y,z);
         if (num_usan>=1) center_val1=usan_vol1.value(x,y,z);
         if (num_usan>=2) center_val2=usan_vol2.value(x,y,z);
	 for(int mz=MAX(zmin,lz); mz<=zmax; mz++) 
	   for(int my=MAX(ymin,ly); my<=ymax; my++) 
	     for(int mx=MAX(xmin,lx); mx<=xmax; mx++) 
	       if ((factor=(float)kernel.value(mx-xmin,my-ymin,mz-zmin)))
	       { 
		 if (num_usan==0) factor*= lut[0][(int)((source.value(mx,my,mz)-center_val)/range)];
		 else factor*= lut[1][(int)((usan_vol1.value(mx,my,mz)-center_val1)/range1)]; 
		 if (num_usan>=2) factor*=lut[2][(int)((usan_vol2.value(mx,my,mz)-center_val2)/range2)]; 
		 num+=source.value(mx,my,mz) * factor;
		 denom+=factor;
	       }
	 if (num_usan>=1) usan_area->value(x,y,z)=(T) denom;
	     if (use_median && denom<1.5)
             {	  
               int count=1;
               for(int x2=MAX(x-1,lx);x2<=MIN(x+1,ux);x2++)
		 for(int y2=MAX(y-1,ly);y2<=MIN(y+1,uy);y2++)
		   for(int z2=MAX(z-medoffst,lz);z2<=MIN(z+medoffst,uz);z2++)
		     if ( (x2-x) || (y2-y) || (z2-z) ) mediankernel(count++)=source.value(x2,y2,z2);
	       ColumnVector subkernel = mediankernel.SubMatrix(1,count-1,1,1);
	       SortAscending(subkernel);        
	       result(x,y,z) = (T)((subkernel(count/2)+subkernel((count+1)/2))/2.0);    
	     }
             else result.value(x,y,z)=(T) (num/denom);
      }
  for(int i=0;i<=num_usan;i++) delete[] (lut[i]-lutsize);
  delete[] lut;
  return result;
}

  template <class T, class S>
  volume<T> convolve(const volume<T>& source, const volume<S>& kernel)
    { 
      extrapolation oldex = source.getextrapolationmethod();
      int offset=0;
      if ((oldex==boundsassert) || (oldex==boundsexception)) 
	{ source.setextrapolationmethod(constpad); }
      volume<T> result(source);
      if (    (( (kernel.maxz() - kernel.minz()) % 2)==1) || 
	      (( (kernel.maxy() - kernel.miny()) % 2)==1) ||
	      (( (kernel.maxx() - kernel.minx()) % 2)==1) ) 
	{
	  cerr << "WARNING:: Off-centre convolution being performed as kernel"
	       << " has even dimensions" << endl;
          //offset=2;
	  //offset gives correction to convolve to match results with fft for even kernel
          //not even kernel with normalise (e.g. -fmean in fslmaths) still has edge problems
	}
      int midx, midy, midz;
      midz=(kernel.maxz() - kernel.minz())/2 + offset;
      midy=(kernel.maxy() - kernel.miny())/2 + offset;
      midx=(kernel.maxx() - kernel.minx())/2 + offset;

      float val;
      for (int z=result.minz(); z<=result.maxz(); z++) {
	for (int y=result.miny(); y<=result.maxy(); y++) {
	  for (int x=result.minx(); x<=result.maxx(); x++) {
	    val=0.0;
	    for (int mz=kernel.minz(); mz<=kernel.maxz(); mz++) {
	      for (int my=kernel.miny(); my<=kernel.maxy(); my++) {
		for (int mx=kernel.minx(); mx<=kernel.maxx(); mx++) {
		  val+=source(x+mx-midx,y+my-midy,z+mz-midz) * kernel(mx,my,mz);
		}
	      }
	    }
	    result(x,y,z)=(T) val;
	  }
	}
      }
      source.setextrapolationmethod(oldex);
      return result;
    }

  template <class T, class S, class M>
  volume<T> convolve(const volume<T>& source, const volume<S>& kernel, 
		     const volume<M>& mask, bool ignoremask, bool renormalise)
    {
      if (!ignoremask && !samesize(mask, source))  
	imthrow("convolve: mask and source are not the same size",10);
      
      if (    (( (kernel.maxz() - kernel.minz()) % 2)==1) || 
	      (( (kernel.maxy() - kernel.miny()) % 2)==1) ||
	      (( (kernel.maxx() - kernel.minx()) % 2)==1) ) 
	{
	  cerr << "WARNING:: Off-centre convolution being performed as kernel"
	       << " has even dimensions" << endl; 
	}
      volume<T> result(source);
      int midx, midy, midz;
      int lz,uz,lx,ux,ly,uy;
      lz=result.minz();
      uz=result.maxz();
      lx=result.minx();
      ux=result.maxx();
      ly=result.miny();
      uy=result.maxy();
      midz=(kernel.maxz() - kernel.minz())/2;
      midy=(kernel.maxy() - kernel.miny())/2;
      midx=(kernel.maxx() - kernel.minx())/2;
      for (int z=lz; z<=uz; z++) 
	for (int y=ly; y<=uy; y++) 
	  for (int x=lx; x<=ux; x++)  
	    if (ignoremask || mask(x,y,z)>0.5) 
            {
	      float val(0),norm(0);
              int x3,y3,z3;
              x3=x-midx;
              y3=y-midy;
              z3=z-midz;
	      for (int mz=kernel.minz(); mz<=kernel.maxz(); mz++) 
		for (int my=kernel.miny(); my<=kernel.maxy(); my++) 
		  for (int mx=kernel.minx(); mx<=kernel.maxx(); mx++) 
                  {
                    int x2,y2,z2;
                    x2=x3+mx;
                    y2=y3+my;
                    z2=z3+mz;
 		    if ((ignoremask && (x2<=ux && x2>=lx && y2<=uy && y2>=ly && z2<=uz && z2>=lz))  || (!ignoremask && mask(x2,y2,z2)>0.5)) 
                    {
		      val+=source.value(x2,y2,z2) * kernel.value(mx,my,mz);
		      norm+=kernel.value(mx,my,mz);
		    }
		  }
	      if (renormalise && fabs(norm)>1e-12) result.value(x,y,z)=(T) (val/norm);
	      else result.value(x,y,z)=(T) val; 
	    }
      return result;
    }

  template <class T>
  volume<T> convolve_separable(const volume<T>& source, 
			       const ColumnVector& kernelx, 
			       const ColumnVector& kernely,
			       const ColumnVector& kernelz)
    {
      volume<T> result(source);
      volume<double> kerx(kernelx.Nrows(),1,1);
      volume<double> kery(1,kernely.Nrows(),1);
      volume<double> kerz(1,1,kernelz.Nrows());
      for (int n=1; n<=kernelx.Nrows(); n++)  kerx.value(n-1,0,0) = kernelx(n);
      for (int n=1; n<=kernely.Nrows(); n++)  kery.value(0,n-1,0) = kernely(n);
      for (int n=1; n<=kernelz.Nrows(); n++)  kerz.value(0,0,n-1) = kernelz(n);
      result = convolve(result,kerx);
      result = convolve(result,kery);
      result = convolve(result,kerz);
      return result;
    }

  template <class T, class M>
  volume<T> convolve_separable(const volume<T>& source, 
			       const ColumnVector& kernelx, 
			       const ColumnVector& kernely,
			       const ColumnVector& kernelz,
			       const volume<M>& mask, bool ignoremask, bool renormalise)
    {
      volume<T> result(source);
      volume<double> kerx(kernelx.Nrows(),1,1);
      volume<double> kery(1,kernely.Nrows(),1);
      volume<double> kerz(1,1,kernelz.Nrows());
      for (int n=1; n<=kernelx.Nrows(); n++)  kerx.value(n-1,0,0) = kernelx(n);
      for (int n=1; n<=kernely.Nrows(); n++)  kery.value(0,n-1,0) = kernely(n);
      for (int n=1; n<=kernelz.Nrows(); n++)  kerz.value(0,0,n-1) = kernelz(n);
      result = convolve(result,kerx,mask,ignoremask,renormalise);
      result = convolve(result,kery,mask,ignoremask,renormalise);
      result = convolve(result,kerz,mask,ignoremask,renormalise);
      return result;
    }

 ///////////////////////////////////////////////////////////////////////////
  // GENERAL DILATION INCLUDING MEDIAN FILTERING

  template <class T, class S>
  volume<T> morphfilter(const volume<T>& source, const volume<S>& kernel,
			const string& filtertype)
    {
      extrapolation oldex = source.getextrapolationmethod();
      if ((oldex==boundsassert) || (oldex==boundsexception)) 
	{ source.setextrapolationmethod(constpad); }
      volume<T> result(source);
      result = 0;

      int nker;
      {
	volume<S> dummy(kernel);
	dummy.binarise((S)0.5);  //new cast to avoid warnings for int-type templates when compiling
	nker = (int) dummy.sum();
      }
      
      int midx, midy, midz;
      if (    (( (kernel.maxz() - kernel.minz()) % 2)==1) || 
	      (( (kernel.maxy() - kernel.miny()) % 2)==1) ||
	      (( (kernel.maxx() - kernel.minx()) % 2)==1) ) 
	{
	  cerr << "WARNING:: Off-centre morphfilter being performed as kernel"
	       << " has even dimensions" << endl;
	}
      midz=(kernel.maxz() - kernel.minz())/2;
      midy=(kernel.maxy() - kernel.miny())/2;
      midx=(kernel.maxx() - kernel.minx())/2;
      int count=1;

      for (int z=result.minz(); z<=result.maxz(); z++) {
	for (int y=result.miny(); y<=result.maxy(); y++) {
	  for (int x=result.minx(); x<=result.maxx(); x++) {
	    ColumnVector vals(nker);
	    count=1;
	    for (int mz=Max(kernel.minz(),result.minz()-z+midz); 
		 mz<=Min(kernel.maxz(),result.maxz()-z+midz); mz++) {
	      for (int my=Max(kernel.miny(),result.miny()-y+midy); 
		   my<=Min(kernel.maxy(),result.maxy()-y+midy); my++) {
		for (int mx=Max(kernel.minx(),result.minx()-x+midx); 
		     mx<=Min(kernel.maxx(),result.maxx()-x+midx); mx++) {
		  if (kernel(mx,my,mz)>0.5) {
		    if ((filtertype!="dilateM" && filtertype!="dilateD") || source(x+mx-midx,y+my-midy,z+mz-midz)) vals(count++) = source(x+mx-midx,y+my-midy,z+mz-midz);
		  }
		}
	      }
	    }
	    if (count>1) {
	      ColumnVector littlevals;
	      littlevals = vals.SubMatrix(1,count-1,1,1);
	      if (filtertype=="median") {  
		SortAscending(littlevals);                         //count/2 works for odd kernel (even count) count+1 gives edge compatibility
		result(x,y,z) = (T)littlevals(Max(1,(count+1)/2)); //with steves IP, otherwise gives the IP median-1 element 
	      } else if ((filtertype=="max") || (filtertype=="dilate")) result(x,y,z) = (T)littlevals.Maximum();
	        else if ((filtertype=="min") || (filtertype=="erode"))   result(x,y,z) = (T)littlevals.Minimum();
                else if (filtertype=="erodeS") {if (source(x,y,z)!=0 && littlevals.Minimum()==0) result(x,y,z) = 0; else result(x,y,z) = source(x,y,z);}
                else if (filtertype=="dilateM") {if (source(x,y,z)==0) result(x,y,z) = (T)(littlevals.Sum()/--count); else result(x,y,z) = source(x,y,z);}
                else if (filtertype=="dilateD") {if (source(x,y,z)==0){
		SortDescending(littlevals); 
                double max=littlevals(1);
                int maxn=1;
                double current=littlevals(1);
                int currentn=1;
                for(int i=2;i<count;i++)
		{
		  if (littlevals(i)==current) currentn++;
                  else
		  {
                    current=littlevals(i);
		    if (currentn>maxn)
		    {
                      max=littlevals(i-1);
                      maxn=currentn; 
                    }
                    currentn=1;
                  }
                }
		result(x,y,z) = (T)max;} else result(x,y,z) = source(x,y,z);
		}
	        else imthrow("morphfilter:: Filter type " + filtertype + "unsupported",7);
	    }   else result(x,y,z) = source(x,y,z);  // THE DEFAULT CASE (leave alone)
	  }
	}
      }
      source.setextrapolationmethod(oldex);
      return result;
    }



 ///////////////////////////////////////////////////////////////////////////
  // RESAMPLE

  template <class T>
  volume4D<T> subsample_by_2(const volume4D<T>& refvol, bool centred)
  {
      volume4D<T> temp_volume;
      for (int t=0; t<refvol.tsize(); t++) {
	temp_volume.addvolume(subsample_by_2(refvol[t],centred));
      }
      // cannot copy the properties as it resets voxel size - but is anything missing
      //     from the 4D properties?
      return temp_volume;
  }

  ///////////////////////////////////////////////////////////////////////////
  // BLURRING

  template <class T>
  volume<T> blur(const volume<T>& source, const ColumnVector& resel_size)
    {
      ColumnVector bmaskx, bmasky, bmaskz;
      make_blur_mask(bmaskx,resel_size(1),source.xdim());
      make_blur_mask(bmasky,resel_size(2),source.ydim());
      make_blur_mask(bmaskz,resel_size(3),source.zdim());
      return convolve_separable(source,bmaskx,bmasky,bmaskz);
    }


  template <class T>
  volume<T> blur(const volume<T>& source, float iso_resel_size)
    {
      ColumnVector resel_size(3);
      resel_size = iso_resel_size;
      return blur(source,resel_size);
    }

  template <class T>
  volume4D<T> blur(const volume4D<T>& source, const ColumnVector& resel_size)
  {
    volume4D<T> dest(source);
    for (int t=source.mint(); t<=source.maxt(); t++) {
      dest[t] = blur(source[t],resel_size);
    }
    return dest;
  }

  template <class T>
  volume4D<T> blur(const volume4D<T>& source, float iso_resel_size)
  {
    volume4D<T> dest(source);
    for (int t=source.mint(); t<=source.maxt(); t++) {
      dest[t] = blur(source[t],iso_resel_size);
    }
    return dest;
  }



  template <class T>
  volume<T> smooth(const volume<T>& source, float sigma_mm)
    {
      float sigmax, sigmay, sigmaz;
      sigmax = sigma_mm/source.xdim();
      sigmay = sigma_mm/source.ydim();
      sigmaz = sigma_mm/source.zdim();
      int nx=((int) (sigmax-0.001))*2 + 3;
      int ny=((int) (sigmay-0.001))*2 + 3;
      int nz=((int) (sigmaz-0.001))*2 + 3;
      ColumnVector kernelx, kernely, kernelz;
      kernelx = gaussian_kernel1D(sigmax,nx);
      kernely = gaussian_kernel1D(sigmay,ny);
      kernelz = gaussian_kernel1D(sigmaz,nz);
      return convolve_separable(source,kernelx,kernely,kernelz);
    }




  template <class T>
  volume<T> smooth2D(const volume<T>& source, float sigma_mm, int nulldir)
    {
      float sigmax, sigmay, sigmaz;
      sigmax = sigma_mm/source.xdim();
      sigmay = sigma_mm/source.ydim();
      sigmaz = sigma_mm/source.zdim();
      int nx=((int) (sigmax-0.001))*2 + 3;
      int ny=((int) (sigmay-0.001))*2 + 3;
      int nz=((int) (sigmaz-0.001))*2 + 3;
      ColumnVector kernelx, kernely, kernelz, nullker(1);
      kernelx = gaussian_kernel1D(sigmax,nx);
      kernely = gaussian_kernel1D(sigmay,ny);
      kernelz = gaussian_kernel1D(sigmaz,nz);
      nullker = 1;
      if (nulldir==1) {
	return convolve_separable(source,nullker,kernely,kernelz);
      } else if (nulldir==2) {
	return convolve_separable(source,kernelx,nullker,kernelz);
      } else {
	// Smoothing in the x-y plane is the default!
	return convolve_separable(source,kernelx,kernely,nullker);
      }
    }

  template <class T>
  volume4D<T> smooth(const volume4D<T>& source,  float sigma_mm)
  {
    volume4D<T> dest(source);
    for (int t=source.mint(); t<=source.maxt(); t++) {
      dest[t] = smooth(source[t],sigma_mm);
    }
    return dest;
  }


  template <class T>
  volume4D<T> smooth(const volume4D<T>& source,  float sigma_mm, int nulldir)
  {
    volume4D<T> dest(source);
    for (int t=source.mint(); t<=source.maxt(); t++) {
      dest[t] = smooth(source[t],sigma_mm,nulldir);
    }
    return dest;
  }

template <class T, class S>
  int insertpart(volume<T>& v1, const volume<S>& v2)  //N.B. This has superficial similarities
{                                                     //to copyconvert, but is in fact quite  
  for (int z=v2.minz(); z<=v2.maxz(); z++) {          //different...
    for (int y=v2.miny(); y<=v2.maxy(); y++) {
      for (int x=v2.minx(); x<=v2.maxx(); x++) {
	v1(x,y,z)=(T) v2(x,y,z);
      }
    }
  }
  return 0;
}


 template <class T, class S,class U>
volume<S> extractpart(const volume<T>& v1, const volume<S>& v2, const volume<U>& kernel) 
{
  volume<S> vout=v2;
  vout = (S) 0.0;
  int kxoff = (kernel.xsize()-1)/2;
  int kyoff = (kernel.ysize()-1)/2;
  int kzoff = (kernel.zsize()-1)/2;
  for (int z=v2.minz(); z<=v2.maxz(); z++) {
    for (int y=v2.miny(); y<=v2.maxy(); y++) {
      for (int x=v2.minx(); x<=v2.maxx(); x++) {
	vout(x,y,z)=(S) v1(x+kxoff,y+kyoff,z+kzoff);
      }
    }
  }
  return vout;
}

////////////////////////////////////////////////////////////////////////////
// Efficient FFT-based convolve
  template <class T, class S>
    volume<T> efficient_convolve(const volume<T>& vin, const volume<S>& vker)
{
  bool usefft=true;
    // estimate calculation time for the two methods and pick the best
    float offt = 2 * vin.nvoxels() * fsllog2(2 * vin.nvoxels());
    float osum = (float)vin.nvoxels() * (float)vker.nvoxels();
    // float cast to avoud overflow for large int multiplication
    float scalefactor = 44;  // relative unit operation cost for fft vs sum
    usefft = (osum > offt * scalefactor);
    //cout << usefft << endl;
    if (usefft) {
    int sx = Max(vin.xsize(),vker.xsize())*2;
    int sy = Max(vin.ysize(),vker.ysize())*2;
    int sz = Max(vin.zsize(),vker.zsize())*2;
    complexvolume vif, vkf;
    vif.re().reinitialize(sx,sy,sz);
    //vif.re().copyproperties(vin);     Is this needed check with MJ...
    vif.re() = 0.0;
    vif.im() = vif.re();
    vkf = vif;
    insertpart(vif.re(),vin);
    insertpart(vkf.re(),vker);
    fft3(vif);
    fft3(vkf);
    vif *= vkf;
    ifft3(vif);
    return extractpart(vif.re(),vin,vker);
    } else return convolve(vin,vker);
}

  template <class T, class S> 
    volume4D<T> generic_convolve(const volume4D<T>& source, const volume<S>& kernel, bool seperable, bool renormalise) 
  { 
      volume4D<T> result(source);
      if (seperable)
      {
        volume<double> kerx(kernel.xsize(),1,1);
        volume<double> kery(1,kernel.ysize(),1);
        volume<double> kerz(1,1,kernel.zsize());
        for (int n=0; n<kernel.xsize(); n++)  kerx.value(n,0,0) = kernel(n,kernel.ysize()/2,kernel.zsize()/2);
        for (int n=0; n<kernel.ysize(); n++)  kery.value(0,n,0) = kernel(kernel.xsize()/2,n,kernel.zsize()/2);
        for (int n=0; n<kernel.zsize(); n++)  kerz.value(0,0,n) = kernel(kernel.xsize()/2,kernel.ysize()/2,n);
        volume<T> mask(1,1,1);
        for (int t=source.mint(); t<=source.maxt(); t++)  result[t] = convolve(result[t],kerx,mask,true,renormalise);
        for (int t=source.mint(); t<=source.maxt(); t++)  result[t] = convolve(result[t],kery,mask,true,renormalise);
        for (int t=source.mint(); t<=source.maxt(); t++)  result[t] = convolve(result[t],kerz,mask,true,renormalise);
      }
      else 
      {
        volume<S> norm_kernel(kernel);
        if (kernel.sum()) norm_kernel/=kernel.sum();
	for (int t=source.mint(); t<=source.maxt(); t++) result[t]=efficient_convolve(source[t],norm_kernel);
        result.copyproperties(source);
        if(renormalise)
	{
	  volume4D<T> unitary_mask(source);
          unitary_mask=1;
          for (int t=source.mint(); t<=source.maxt(); t++) unitary_mask[t]=efficient_convolve(unitary_mask[t],norm_kernel);
          result/=unitary_mask;
        }
      }
    return result; 
  } 


  ///////////////////////////////////////////////////////////////////////////
  // GRADIENT


  template <class T>
  volume<float> gradient(const volume<T>& source)
    {
      volume<float> maskx,masky,maskz;
      make_grad_masks(maskx,masky,maskz);
      volume<float> grad(source);
      float valx, valy, valz;
      int midx, midy, midz;
      midz=maskx.xsize()/2;
      midy=maskx.ysize()/2;
      midx=maskx.zsize()/2;
      for (int z=0; z<grad.zsize(); z++) {
	for (int y=0; y<grad.ysize(); y++) {
	  for (int x=0; x<grad.xsize(); x++) {
	    valx=0.0; valy=0.0; valz=0.0;
	    for (int mz=-midz; mz<=midz; mz++) {
	      for (int my=-midy; my<=midy; my++) {
		for (int mx=-midx; mx<=midx; mx++) {
		  valx+=source(x+mx,y+my,z+mz) * maskx(mx+midx,my+midy,mz+midz);
		  valy+=source(x+mx,y+my,z+mz) * masky(mx+midx,my+midy,mz+midz);
		  valz+=source(x+mx,y+my,z+mz) * maskz(mx+midx,my+midy,mz+midz);
		}
	      }
	    }
	    grad(x,y,z)=sqrt(Sqr(valx) + Sqr(valy) + Sqr(valz));
	  }
	}
      }
      return grad;
    }



   template <class T>
   void gradient(const volume<T>& source,volume4D<float>& grad)
    {
      volume<float> maskx,masky,maskz;
      make_grad_masks(maskx,masky,maskz);
      grad.reinitialize(source.xsize(),source.ysize(),source.zsize(),3);
      copybasicproperties(source,grad[0]);
      float valx, valy, valz;
      int midx, midy, midz;
      midz=maskx.xsize()/2;
      midy=maskx.ysize()/2;
      midx=maskx.zsize()/2;
      for (int z=0; z<grad.zsize(); z++) {
	for (int y=0; y<grad.ysize(); y++) {
	  for (int x=0; x<grad.xsize(); x++) {
	    valx=0.0; valy=0.0; valz=0.0;
	    for (int mz=-midz; mz<=midz; mz++) {
	      for (int my=-midy; my<=midy; my++) {
		for (int mx=-midx; mx<=midx; mx++) {
		  valx+=source(x+mx,y+my,z+mz) * maskx(mx+midx,my+midy,mz+midz);
		  valy+=source(x+mx,y+my,z+mz) * masky(mx+midx,my+midy,mz+midz);
		  valz+=source(x+mx,y+my,z+mz) * maskz(mx+midx,my+midy,mz+midz);
		}
	      }
	    }
	    grad(x,y,z,0)=valx;
	    grad(x,y,z,1)=valy;
	    grad(x,y,z,2)=valz;
	  }
	}
      }
      
    }



  template <class T>
  volume4D<float> lrxgrad(const volume<float>& im, const volume<T>& mask) 
  {
    // calculates separate left and right gradients (with copying when at
    //  borders of the mask) or zero for points outside of the mask
    // returns the two as part of a volume4D (vol[0] = left grad, vol[1] = right grad)
    volume4D<float> grad;
    if (!samesize(im,mask)) imthrow("Mask and image not the same size",20);
    grad.addvolume(im);
    grad.addvolume(im);
    for (int z=im.minz(); z<=im.maxz(); z++) {
      for (int y=im.miny(); y<=im.maxy(); y++) {
	for (int x=im.minx(); x<=im.maxx(); x++) {
	  // defaults, only overridden if inside mask and can calculate the gradients
	  grad[0](x,y,z) = 0;
	  grad[1](x,y,z) = 0;
	  if (mask(x,y,z)>0.5) {
	    // left gradient
	    if (x>im.minx()) {
	      if (mask(x-1,y,z)>0.5) {
		grad[0](x,y,z) = im(x,y,z) - im(x-1,y,z);
	      }
	    }
	    // right gradient
	    if (x<im.maxx()) {
	      if (mask(x+1,y,z)>0.5) {
		grad[1](x,y,z) = im(x+1,y,z) - im(x,y,z);
	      } else {
		// if couldn't calculate right grad, then copy left
		grad[1](x,y,z) = grad[0](x,y,z);
	      }
	    }
	    // if couldn't calculate left grad, then copy right
	    if ( (x>im.minx()) && (mask(x-1,y,z)<=0.5) ) {
	      grad[0](x,y,z) = grad[1](x,y,z);
	    }
	    // NB: if couldn't calculate either (but mask>0.5) then it is still 0
	  }
	}
      }
    }
    return grad;
  }


  template <class T>
  volume4D<float> lrygrad(const volume<float>& im, const volume<T>& mask) 
  {
    // calculates separate left and right gradients (with copying when at
    //  borders of the mask) or zero for points outside of the mask
    // returns the two as part of a volume4D (vol[0] = left grad, vol[1] = right grad)
    volume4D<float> grad;
    if (!samesize(im,mask)) imthrow("Mask and image not the same size",20);
    grad.addvolume(im);
    grad.addvolume(im);
    for (int z=im.minz(); z<=im.maxz(); z++) {
      for (int y=im.miny(); y<=im.maxy(); y++) {
	for (int x=im.minx(); x<=im.maxx(); x++) {
	  // defaults, only overridden if inside mask and can calculate the gradients
	  grad[0](x,y,z) = 0;
	  grad[1](x,y,z) = 0;
	  if (mask(x,y,z)>0.5) {
	    // left gradient
	    if (y>im.miny()) {
	      if (mask(x,y-1,z)>0.5) {
		grad[0](x,y,z) = im(x,y,z) - im(x,y-1,z);
	      }
	    }
	    // right gradient
	    if (y<im.maxy()) {
	      if (mask(x,y+1,z)>0.5) {
		grad[1](x,y,z) = im(x,y+1,z) - im(x,y,z);
	      } else {
		// if couldn't calculate right grad, then copy left
		grad[1](x,y,z) = grad[0](x,y,z);
	      }
	    }
	    // if couldn't calculate left grad, then copy right
	    if ( (y>im.miny()) && (mask(x,y-1,z)<=0.5) ) {
	      grad[0](x,y,z) = grad[1](x,y,z);
	    }
	    // NB: if couldn't calculate either (but mask>0.5) then it is still 0
	  }
	}
      }
    }
    return grad;
  }


  template <class T>
  volume4D<float> lrzgrad(const volume<float>& im, const volume<T>& mask) 
  {
    // calculates separate left and right gradients (with copying when at
    //  borders of the mask) or zero for points outside of the mask
    // returns the two as part of a volume4D (vol[0] = left grad, vol[1] = right grad)
    volume4D<float> grad;
    if (!samesize(im,mask)) imthrow("Mask and image not the same size",20);
    grad.addvolume(im);
    grad.addvolume(im);
    for (int z=im.minz(); z<=im.maxz(); z++) {
      for (int y=im.miny(); y<=im.maxy(); y++) {
	for (int x=im.minx(); x<=im.maxx(); x++) {
	  // defaults, only overridden if inside mask and can calculate the gradients
	  grad[0](x,y,z) = 0;
	  grad[1](x,y,z) = 0;
	  if (mask(x,y,z)>0.5) {
	    // left gradient
	    if (z>im.minz()) {
	      if (mask(x,y,z-1)>0.5) {
		grad[0](x,y,z) = im(x,y,z) - im(x,y,z-1);
	      }
	    }
	    // right gradient
	    if (z<im.maxz()) {
	      if (mask(x,y,z+1)>0.5) {
		grad[1](x,y,z) = im(x,y,z+1) - im(x,y,z);
	      } else {
		// if couldn't calculate right grad, then copy left
		grad[1](x,y,z) = grad[0](x,y,z);
	      }
	    }
	    // if couldn't calculate left grad, then copy right
	    if ( (z>im.minz()) && (mask(x,y,z-1)<=0.5) ) {
	      grad[0](x,y,z) = grad[1](x,y,z);
	    }
	    // NB: if couldn't calculate either (but mask>0.5) then it is still 0
	  }
	}
      }
    }
    return grad;
  }


  
  
   template <class T>
     volume4D<T> bandpass_temporal_filter(volume4D<T>& source,double hp_sigma, double lp_sigma)
      { 
	//cout << "hp " << hp_sigma << "lp " << lp_sigma << endl;
         int backwards = 0;
         int hp_mask_size_PLUS, lp_mask_size_PLUS, hp_mask_size_MINUS, lp_mask_size_MINUS;
         double *hp_exp=NULL, *lp_exp=NULL, *array, *array2;
         volume4D<T> result(source);

         if (hp_sigma<=0) hp_mask_size_MINUS=0;
         else hp_mask_size_MINUS=(int)(hp_sigma*3);   /* this isn't a linear filter, so small hard cutoffs at ends don't matter */
         if (!backwards) hp_mask_size_PLUS=hp_mask_size_MINUS;
         else hp_mask_size_PLUS=0;
         if (lp_sigma<=0) lp_mask_size_MINUS=0;
         else lp_mask_size_MINUS=(int)(lp_sigma*5)+2; /* this will be small, so we might as well be careful */
         if (!backwards) lp_mask_size_PLUS=lp_mask_size_MINUS;
         else lp_mask_size_PLUS=0;

         array=new double[source.tsize()+2*lp_mask_size_MINUS];
	 array+=lp_mask_size_MINUS;
         array2=new double[source.tsize()+2*lp_mask_size_MINUS];
         array2+=lp_mask_size_MINUS;

         if (hp_sigma>0)
         {
            hp_exp=new double[hp_mask_size_MINUS+hp_mask_size_PLUS+1];
            hp_exp+=hp_mask_size_MINUS;
            for(int t=-hp_mask_size_MINUS; t<=hp_mask_size_PLUS; t++)
            hp_exp[t] = exp( -0.5 * ((double)(t*t)) / (hp_sigma*hp_sigma) );
         }

         if (lp_sigma>0)
         {
            double total=0;
            lp_exp=new double[lp_mask_size_MINUS+lp_mask_size_PLUS+1];
            lp_exp+=lp_mask_size_MINUS;
            for(int t=-lp_mask_size_MINUS; t<=lp_mask_size_PLUS; t++)
            {
              lp_exp[t] = exp( -0.5 * ((double)(t*t)) / (lp_sigma*lp_sigma) );
              total += lp_exp[t];
            }

            for(int t=-lp_mask_size_MINUS; t<=lp_mask_size_PLUS; t++)
              lp_exp[t] /= total;
         }
         for(int z=0;z<source.zsize();z++)
           for(int y=0;y<source.ysize();y++)	    
	     for(int x=0;x<source.xsize();x++)
	     {
               for(int t=0; t<source.tsize(); t++) array[t] = (double)source.value(x,y,z,t);
               if (hp_sigma>0)
               {
                 int done_c0=0;
                 double c0=0;
                 for(int t=0; t<source.tsize(); t++)
                 {
                    int tt;
                    double c, w, A=0, B=0, C=0, D=0, N=0, tmpdenom;
                    for(tt=MAX(t-hp_mask_size_MINUS,0); tt<=MIN(t+hp_mask_size_PLUS,source.tsize()-1); tt++)
                    {
                      int dt=tt-t;
                      w = hp_exp[dt];
                      A += w * dt;
                      B += w * array[tt];
                      C += w * dt * dt;
                      D += w * dt * array[tt];
                      N += w;
                    }
                    tmpdenom=C*N-A*A;
                    if (tmpdenom!=0)
	            {
	               c = (B*C-A*D) / tmpdenom;
	               if (!done_c0)
	               {
	                 c0=c;
	                 done_c0=1;
	               }
	               array2[t] = c0 + array[t] - c;
	             }
	             else  array2[t] = array[t];
	          }
                  memcpy(array,array2,sizeof(double)*source.tsize());
	       }
	       /* {{{ apply lowpass filter to 1D array */
               if (lp_sigma>0)
               {
          /* {{{ pad array at ends */
                 for(int t=1; t<=lp_mask_size_MINUS; t++)
                 {
                    array[-t]=array[0];
                    array[source.tsize()-1+t]=array[source.tsize()-1];
                 }
                 for(int t=0; t<source.tsize(); t++)
                 { 
                    double total=0;
                    int tt;
                    for(tt=t-lp_mask_size_MINUS; tt<=t+lp_mask_size_PLUS; tt++)   total += array[tt] * lp_exp[tt-t];
                    array2[t] = total;
                 }
                  memcpy(array,array2,sizeof(double)*source.tsize());
	       }
	  /* {{{ write 1D array back to input 4D data */
               for(int t=0; t<source.tsize(); t++) result.value(x,y,z,t)= (T)array[t];
	     }
	 return result;
      }


  ///////////////////////////////////////////////////////////////////////////
  // EDGE DETECTION

  /* detects closed contour edges in a volume as the zero crossings of the 
     Laplacian of a Gaussian (LoG). This is implemented by convolving the 
     data with a kernel formed by subtracting a Gaussian of radius sigma1
     from a second Gaussian kernel of radius sigma2 (sigma1 < sigma2) */

  template <class T>
  volume<T> log_edge_detect(const volume<T>& source, 
			    float sigma1, float sigma2, 
			    float zero_tolerance, bool twodimensional)
    {
      /* zero_tolerance is what we define as an acceptable "zero-crossing" */
      int radius1;
      volume<float> log_kern, temp_kern;
      volume<T> result(source);

      radius1 = (int)(4*sigma2);
      if (twodimensional) {
	log_kern = gaussian_kernel2D(sigma2, radius1);
	temp_kern = gaussian_kernel2D(sigma1, radius1);
      } else {
	log_kern = gaussian_kernel3D(sigma2, radius1);
	temp_kern = gaussian_kernel3D(sigma1, radius1);
      }
      
      log_kern -= temp_kern;

      result = convolve(source, log_kern);
      result.binarise(zero_tolerance);

      return result;
    }



   template <class T>
   volume<T> fixed_edge_detect(const volume<T>& source, float threshold, 
			       bool twodimensional)
   {
     volume<T> result = source;
     int zsize = 3;
     if (twodimensional) zsize=1;
     
     volume<float> log_kern(3,3,zsize);
     log_kern = -1;
     log_kern(1,1,(zsize-1)/2) = 8;

     extrapolation oldex = source.getextrapolationmethod();
     source.setextrapolationmethod(mirror);
     result = convolve(source, log_kern);
     source.setextrapolationmethod(oldex);
     result.binarise(threshold);

     return result;
   }



  ///////////////////////////////////////////////////////////////////////////
  // EDGE STRENGTHEN (from avwmaths)
  template <class T>
  volume4D<T> edge_strengthen(const volume4D<T>& source)
  {
    float tmpf=2*sqrt(1/(pow((double)source.xdim(),2.0)) + 1/(pow((double)source.ydim(),2.0)) + 1/(pow((double)source.zdim(),2.0)));
       volume4D<T> result;
       result=source;
       if (source.zsize()>2)
       {
          for(int t=0;t<source.tsize();t++)           
           for(int z=1;z<source.zsize()-1;z++)
             for(int y=1;y<source.ysize()-1;y++)	    
	       for(int x=1;x<source.xsize()-1;x++)
	       {
                 double temp1 = pow(double(source.value(x,y,z+1,t)-source.value(x,y,z-1,t)),2.0)/ pow((double)source.zdim(),2.0);
                 double temp2 = pow(double(source.value(x,y+1,z,t)-source.value(x,y-1,z,t)),2.0)/ pow((double)source.ydim(),2.0);
		 double temp3 = pow(double(source.value(x+1,y,z,t)-source.value(x-1,y,z,t)),2.0)/ pow((double)source.xdim(),2.0);
		 result(x,y,z,t)=(T)(sqrt(temp1+temp2+temp3)/tmpf);
	      }
       } 
       else 
       {
         for(int t=0;t<source.tsize();t++)           
           for(int z=0;z<source.zsize();z++)
             for(int y=1;y<source.ysize()-1;y++)	    
	       for(int x=1;x<source.xsize()-1;x++)
	       {
                 double temp1=pow(double(source.value(x,y+1,z,t)-source.value(x,y-1,z,t)),2.0)/ pow((double)source.ydim(),2.0);            
		 double temp2=pow(double(source.value(x+1,y,z,t)-source.value(x-1,y,z,t)),2.0)/ pow((double)source.xdim(),2.0);
		 result(x,y,z,t) = (T)(sqrt (temp1+temp2) / tmpf);
               }
       }
       return result;
  }


  ///////////////////////////////////////////////////////////////////////////
  // CONNECTED COMPONENTS

  // support functions for connected components
  
  int find_first_nonzero(const Matrix& mat);
 
  void addpair2set(int x, int y, std::vector<int>& sx, std::vector<int>& sy);

 void relabel_components_uniquely(volume<int>& labelvol, 
				   const std::vector<int>& equivlista,
				   const std::vector<int>& equivlistb, ColumnVector& clustersizes); 

  void relabel_components_uniquely(volume<int>& labelvol, 
				   const std::vector<int>& equivlista,
				   const std::vector<int>& equivlistb); 

  ////////////////////////////////////////////////////////////////////////////

  template <class T>
  void nonunique_component_labels(const volume<T>& vol, 
				  volume<int>& labelvol, 
				  std::vector<int>& equivlista,
				  std::vector<int>& equivlistb,
				  int numconnected)
    {
      copyconvert(vol,labelvol);
      labelvol = 0;
  
      int labelnum=1;

      equivlista.erase(equivlista.begin(),equivlista.end());
      equivlistb.erase(equivlistb.begin(),equivlistb.end());

      for (int z=vol.minz(); z<=vol.maxz(); z++) {
	for (int y=vol.miny(); y<=vol.maxy(); y++) {
	  for (int x=vol.minx(); x<=vol.maxx(); x++) {
	    T val = vol(x,y,z);
	    if (val>0.5) {  // The eligibility test
	      int lval = labelvol(x,y,z);
	      for (int znew=z-1; znew<=z; znew++) {
		int ystart=y, yend=y;
		if (numconnected==6) {
		  ystart += -1-znew+z;
		} else {
		  ystart += -1;
		  yend += z-znew;
		}
		for (int ynew=ystart; ynew<=yend; ynew++) {
		  int xstart=x, xend=x;
		  if (numconnected==6) {
		    xstart += -1 -(znew-z) -(ynew-y) -(znew-z)*(ynew-y);
		    xend = xstart;
		  } else if (numconnected==18) {
		    xstart += -1 -(znew-z)*std::abs(ynew-y);
		    xend = xstart + 2*(znew - z + std::abs(ynew-y));
		  } else {
		    xstart += -1;
		    xend += 1 -2*(znew-z+1)*(ynew-y+1);
		  }
		  for (int xnew=xstart; xnew<=xend; xnew++) {
		    if ( (xnew>=vol.minx()) && (ynew>=vol.miny()) 
			 && (znew>=vol.minz())
			 && (MISCMATHS::round(vol(xnew,ynew,znew)-val)==0) ) { 
		      // Binary relation
		      int lval2 = labelvol(xnew,ynew,znew);
		      if (lval != lval2) {
			if (lval!=0) {
			  addpair2set(lval2,lval,equivlista,equivlistb);
			}
			labelvol(x,y,z) = lval2;
			lval = lval2;
		      }
		    }
		  }
		}
	      }
	      if (lval==0) {
		labelvol(x,y,z) = labelnum;
		labelnum++;
	      }
	    }
	  }
	}
      }
      
    }


  template <class T>
  void nonunique_component_labels(const volume<T>& vol, 
				  const volume<T>& mask,
				  volume<int>& labelvol, 
				  std::vector<int>& equivlista,
				  std::vector<int>& equivlistb,
				  bool (*binaryrelation)(T , T),
				  int numconnected)
    {
      copyconvert(vol,labelvol);
      labelvol = 0;
  
      int labelnum=1;

      equivlista.erase(equivlista.begin(),equivlista.end());
      equivlistb.erase(equivlistb.begin(),equivlistb.end());

      for (int z=vol.minz(); z<=vol.maxz(); z++) {
	for (int y=vol.miny(); y<=vol.maxy(); y++) {
	  for (int x=vol.minx(); x<=vol.maxx(); x++) {
	    if (mask(x,y,z)>0.5) {  // The eligibility test
	      int lval = labelvol(x,y,z);
	      int val = vol(x,y,z);
	      for (int znew=z-1; znew<=z; znew++) {
		int ystart=y, yend=y;
		if (numconnected==6) {
		  ystart += -1-znew+z;
		} else {
		  ystart += -1;
		  yend += z-znew;
		}
		for (int ynew=ystart; ynew<=yend; ynew++) {
		  int xstart=x, xend=x;
		  if (numconnected==6) {
		    xstart += -1 -(znew-z) -(ynew-y) -(znew-z)*(ynew-y);
		    xend = xstart;
		  } else if (numconnected==18) {
		    xstart += -1 -(znew-z)*std::abs(ynew-y);
		    xend = xstart + 2*(znew - z + std::abs(ynew-y));
		  } else {
		    xstart += -1;
		    xend += 1 -2*(znew-z+1)*(ynew-y+1);
		  }
		  for (int xnew=xstart; xnew<=xend; xnew++) {
		    if ( (xnew>=vol.minx()) && (ynew>=vol.miny()) 
			 && (znew>=vol.minz())
			 && ((*binaryrelation)(vol(xnew,ynew,znew),val)) ) { 
		      // Binary relation
		      int lval2 = labelvol(xnew,ynew,znew);
		      if (lval != lval2) {
			if (lval!=0) {
			  addpair2set(lval2,lval,equivlista,equivlistb);
			}
			labelvol(x,y,z) = lval2;
			lval = lval2;
		      }
		    }
		  }
		}
	      }
	      if (lval==0) {
		labelvol(x,y,z) = labelnum;
		labelnum++;
	      }
	    }
	  }
	}
      }
      
    }


  template <class T>
  volume<int> connected_components(const volume<T>& vol, ColumnVector& clustersize, int numconnected)
    {
      volume<int> labelvol;
      copyconvert(vol,labelvol);
      std::vector<int> equivlista, equivlistb;
      nonunique_component_labels(vol,labelvol,
				 equivlista,equivlistb,numconnected);
      relabel_components_uniquely(labelvol,equivlista,equivlistb,clustersize);
      return labelvol;
    }


  template <class T>
  volume<int> connected_components(const volume<T>& vol, 
                                   const volume<T>& mask, 
                                   bool (*binaryrelation)(T , T), ColumnVector& clustersize)
    {
      volume<int> labelvol;
      copyconvert(vol,labelvol);
      std::vector<int> equivlista, equivlistb;
      nonunique_component_labels(vol,mask,labelvol,equivlista,equivlistb,
				 binaryrelation);
      relabel_components_uniquely(labelvol,equivlista,equivlistb,clustersize);
      return labelvol;
    }
  
  
  template <class T>
  volume<int> connected_components(const volume<T>& vol, int numconnected){
    ColumnVector clustersize;
    return connected_components(vol,clustersize,numconnected);
  }

  
  template <class T>
  volume<int> connected_components(const volume<T>& vol,const volume<T>& mask, 
                                   bool (*binaryrelation)(T , T)){
    ColumnVector clustersize;
    return connected_components(vol,mask,binaryrelation,clustersize);
  }


//////////////////////////////////// distancemap-related functions //////////////////////


class rowentry { public: int x; int y; int z; float d; } ;

bool rowentry_lessthan(const rowentry& r1, const rowentry& r2);

template <class T>
class distancemapper {
private:
  const volume<T> &bvol;
  const volume<T> &mask;
  vector<rowentry> schedule;
  Matrix octantsign;
public:
  distancemapper(const volume<T>& binaryvol, const volume<T>& maskvol);
  ~distancemapper();
  volume<float> distancemap();
  volume4D<float> sparseinterpolate(const volume4D<float>& values, 
				    const string& interpmethod="general");
private:
  int setup_globals();  
  int find_nearest(int x, int y, int z, int& x1, int& y1, int& z1, 
		   bool findav, ColumnVector& localav, const volume4D<float>& vals);
  int find_nearest(int x, int y, int z, int& x1, int& y1, int& z1);
  int create_distancemap(volume4D<float>& vout, const volume4D<float>& valim,
			  const string& interpmethod="none");
};

template <class T>
distancemapper<T>::distancemapper(const volume<T>& binaryvol, const volume<T>& maskvol) :
  bvol(binaryvol), mask(maskvol)
{
  if (!samesize(bvol,mask)) imthrow("Mask and image not the same size",20);
  octantsign.ReSize(8,3);
  setup_globals();
}

template <class T>
distancemapper<T>::~distancemapper()
{
  // destructors for the private members should do the job
}

template <class T>
int distancemapper<T>::setup_globals()
{
  // octantsign gives the 8 different octant sign combinations for coord
  //  offsets
  int row=1;
  for (int p=-1; p<=1; p+=2) {
    for (int q=-1; q<=1; q+=2) {
      for (int r=-1; r<=1; r+=2) {
	octantsign(row,1)=p;
	octantsign(row,2)=q;
	octantsign(row,3)=r;
	row++;
      }
    }
  }

  // construct list of displacements (in one octant) in ascending
  // order of distance
  for (int z=bvol.minz(); z<=bvol.maxz(); z++) {
    for (int y=bvol.miny(); y<=bvol.maxy(); y++) {
      for (int x=bvol.minx(); x<=bvol.maxx(); x++) {
	if ( (x==0) && (y==0) && (z==0) ) 
         { // do nothing 
         } else {
	   rowentry newrow;
	   newrow.x=x;
	   newrow.y=y;
	   newrow.z=z;
	   float d2 = norm2sq(x*bvol.xdim(),y*bvol.ydim(),z*bvol.zdim());
	   newrow.d=d2;
	   schedule.push_back(newrow);
         }
      }
    }
  }

  // sort schedule to get ascending d2
  sort(schedule.begin(),schedule.end(),NEWIMAGE::rowentry_lessthan);
  return 0;
}

// findav determines whether to do interpolation calculations or just
//  return the location only
template <class T>
int distancemapper<T>::find_nearest(int x, int y, int z, int& x1, int& y1, int& z1, 
				 bool findav, ColumnVector& localav,
				 const volume4D<float>& vals)
{
  float sumw=0.0, mindist=0.0, maxdist=0.0, weight;
  ColumnVector sumvw;
  if (findav) { 
    localav.ReSize(vals.tsize()); 
    localav=0.0;
    sumvw.ReSize(vals.tsize());
    sumvw=0.0;
  }
  for (int r=0; r<(int) schedule.size(); r++) 
  {
    for (int p=1; p<=8; p++) 
      {
	int dx=MISCMATHS::round(schedule[r].x*octantsign(p,1));
	int dy=MISCMATHS::round(schedule[r].y*octantsign(p,2));
	int dz=MISCMATHS::round(schedule[r].z*octantsign(p,3));
	if (bvol.in_bounds(x+dx,y+dy,z+dz) && (bvol.value(x+dx,y+dy,z+dz)>0.5)) 
	  { 
	    x1=x+dx; y1=y+dy; z1=z+dz; 
	    if (!findav) { 
	      return 0;
	    } else {
	      if (mindist==0.0) {  // first time a point is encountered
		mindist=schedule[r].d;
		// select distance band to average over (farther -> more)
		maxdist=MISCMATHS::Max(mindist+1.0,mindist*1.5);
	      } else {
		if ( (schedule[r].d>maxdist) )
		  {  // stop after maxdist reached
		    localav=sumvw/sumw;
		    return 0;
		  }
	      }
	      weight = mindist/schedule[r].d;
	      sumw += weight;
	      for (int t=0; t<vals.tsize(); t++) {
		sumvw(t+1) += weight * vals.value(x+dx,y+dy,z+dz,t);
	      }
	    }
	  }
      }
  }
  // return furtherest point (in error)
  x1= x + MISCMATHS::round(schedule.back().x);
  y1= y + MISCMATHS::round(schedule.back().y);
  z1= z + MISCMATHS::round(schedule.back().z);
  if (!findav) { if (sumw>0) { localav=sumvw/sumw; return 0; } else localav=0.0; }
  return 1;
}

template <class T>
int distancemapper<T>::find_nearest(int x, int y, int z, int& x1, int& y1, int& z1)
{ 
  ColumnVector dummy;
  volume4D<float> dummyvol;
  return this->find_nearest(x,y,z,x1,y1,z1,false,dummy,dummyvol);
}


// create the distance map as vout
// if return_distance is false then interpolate the value of the input volume
//  at the output location rather than store the distance to it
template <class T>
int distancemapper<T>::create_distancemap(volume4D<float>& vout, const volume4D<float>& valim,
					  const string& interpmethod)
{
  int x1, y1, z1;
  ColumnVector localav;
  int interp=0;
  if ((interpmethod=="nn") || (interpmethod=="nearestneighbour")) interp=1;
  if (interpmethod=="general") interp=2;
  if (interp>0) { vout = valim; } else { vout = bvol; vout *= 0.0f; }
  if ((interp>0) && (!samesize(bvol,valim[0])))
    { imthrow("Binary image and interpolant not the same size",21); }
  for (int z=vout.minz(); z<=vout.maxz(); z++) {
    for (int y=vout.miny(); y<=vout.maxy(); y++) {
      for (int x=vout.minx(); x<=vout.maxx(); x++) {
	if (mask(x,y,z)>((T) 0.5)) {
	  if (interp>=2) {
	    find_nearest(x,y,z,x1,y1,z1,true,localav,valim);
	  } else {
	    find_nearest(x,y,z,x1,y1,z1);
	  }
	  switch (interp) {
	  case 2:
	    for (int t=0;t<valim.tsize();t++) { vout(x,y,z,t)=localav(t+1); }
	    break;
	  case 1:
	    for (int t=0;t<valim.tsize();t++) { vout(x,y,z,t)=valim(x1,y1,z1,t); }
	    break;
	  case 0:
	  default:
	    vout(x,y,z,0)=sqrt(norm2sq((x1-x)*bvol.xdim(),
				       (y1-y)*bvol.ydim(),(z1-z)*bvol.zdim()));
	  }
	}
      }
    }
  }
  return 0;
}


template <class T>
volume<float> distancemapper<T>::distancemap()
{
  volume4D<float> dmap;
  create_distancemap(dmap,dmap,"none");
  return dmap[0];
}

template <class T>
volume4D<float> distancemapper<T>::sparseinterpolate(const volume4D<float>& values, 
						     const string& interpmethod)
{
  volume4D<float> vout;
  create_distancemap(vout,values,interpmethod);
  return vout;
}


template <class T>
volume<float> distancemap(const volume<T>& binaryvol)
{
  volume<T> mask;
  mask = ((T) 1) - binarise(binaryvol,((T) 0.5));
  return distancemap(binaryvol,mask);
}


template <class T>
volume<float> distancemap(const volume<T>& binaryvol, const volume<T>& mask)
{
  distancemapper<T> dmapper(binaryvol,mask);
  return dmapper.distancemap();
}

template <class T>
volume4D<float> sparseinterpolate(const volume4D<T>& sparsesamps, const volume<T>& mask,
				  const string& interpmethod)
{
  // can have "general" or "nearestneighbour" (or "nn") for interpmethod
  volume<T> invmask;
  invmask=((T) 1) - mask;
  distancemapper<T> dmapper(mask,invmask);
  return dmapper.sparseinterpolate(sparsesamps,interpmethod);
}



//////////////////////////////////// TFCE-related functions //////////////////////

template <class T>
void tfce_orig_slow(volume<T>& VolIntn, float H, float E, int NumConn, float minT, float deltaT)
{
  float maxval=VolIntn.max();
  volume<float> clusterenhance;
  copyconvert(VolIntn,clusterenhance);
  clusterenhance=0;

  if (deltaT==0)
    deltaT = (maxval - minT)/100.0;   // this needs fixing!!

  for (float thresh=minT+deltaT; thresh<=maxval; thresh+=deltaT)
    {
      volume<float> clusters;
      copyconvert(VolIntn,clusters);
      clusters.binarise(thresh);

      ColumnVector clustersizes;  
      volume<int>tmpvol=connected_components(clusters,clustersizes,NumConn);
      clustersizes = pow(clustersizes,E) * pow(thresh,H);
      for(int z=0;z<VolIntn.zsize();z++)
	for(int y=0;y<VolIntn.ysize();y++)	    
	  for(int x=0;x<VolIntn.xsize();x++)
	    if (tmpvol.value(x,y,z)>0)
	      clusterenhance.value(x,y,z) += clustersizes(tmpvol.value(x,y,z));
    }
  copyconvert(clusterenhance,VolIntn);
  return;
}



class VecSort{
 public:
  int Sx, Sy, Sz, Sl;
  double Sv;
  bool operator<(const VecSort& other) const{
    return Sv < other.Sv;
  }
};

//
// minT would normally be 0
// if deltaT is set to 0 it is reset to max/100
//
template <class T>
void tfce(volume<T>& data, float H, float E, int NumConn, float minT, float deltaT)
{
  volume<int> VolLabl; copyconvert(data, VolLabl);
  volume<float> VolEnhn; copyconvert(data, VolEnhn); VolEnhn=0;
  bool doIT=false;
  const int INIT=-1, MASK=-2;
  int curlab=0;
  int pX, pY, pZ, qX, qY, qZ, rX, rY, rZ;
  int FldCntr=0, FldCntri=0, tfceCntr=0, xFC = 0;
  int minX=1, minY=1, minZ=1, maxX=data.maxx()-1, maxY=data.maxy()-1, maxZ=data.maxz()-1;
  int sizeC=maxX*maxY*maxZ;
  int counter=0, edsta[27];
  float maxT=data.max();
  if(deltaT==0) 
    deltaT=maxT/100.0;
  if(deltaT<=0) 
    throw Exception("Error: tfce requires a positive deltaT input.");
  if ( data.xsize() < 3 || data.ysize() < 3 || data.zsize() < 3 )
    throw Exception("Error: tfce currently requires an input with at least 3 voxels extent into each dimension.");
  if( data.max()/deltaT > 10000 )
    cout << "Warning: tfce has detected a large number of integral steps. This operation may require a great deal of time to complete." << endl;
  vector<VecSort> VecSortI(sizeC);
  queue<int> Qx, Qy, Qz;      
  for(int z0=-1; z0<=1; z0++)
    for(int y0=-1; y0<=1; y0++)
      for(int x0=-1; x0<=1; x0++){
	edsta[counter++] = (x0*x0+y0*y0+z0*z0);
	if (edsta[counter-1]<2 && edsta[counter-1]!=0) edsta[counter-1]=6;
	if (edsta[counter-1]<3 && edsta[counter-1]!=0) edsta[counter-1]=18;
	if (edsta[counter-1]<4 && edsta[counter-1]!=0) edsta[counter-1]=26;
	if (edsta[counter-1]==0) edsta[counter-1]=100;
      }
  FldCntr=0;
  for(int z=minZ; z<=maxZ; z++)
    for(int y=minY; y<=maxY; y++)
      for(int x=minX; x<=maxX; x++) {
	float iVal=data.value(x,y,z);
	if( iVal > minT ) {
	  VecSortI[FldCntr].Sx=x; VecSortI[FldCntr].Sy=y; VecSortI[FldCntr].Sz=z;
	  VecSortI[FldCntr++].Sv=iVal;
	}
      }
  sizeC=FldCntr;
  VecSortI.resize(sizeC);
  sort(VecSortI.begin(), VecSortI.end());
  for(float curThr=minT; curThr<(maxT+deltaT);curThr+=deltaT){
    VolLabl = INIT; FldCntr = xFC; curlab = 0;
    while( (VecSortI[FldCntr].Sv<=curThr) && (FldCntr<sizeC) ){
      VecSortI[FldCntr].Sl=0; VecSortI[FldCntr++].Sv=0;
    }
    xFC=FldCntr;
    while( VecSortI[FldCntr].Sv>curThr  && (FldCntr<sizeC) ){
      pX=VecSortI[FldCntr].Sx; pY=VecSortI[FldCntr].Sy; pZ=VecSortI[FldCntr++].Sz;
      VolLabl.value(pX,pY,pZ)=MASK;
    }             
    for(FldCntri=xFC; FldCntri<FldCntr; FldCntri++){
      pX=VecSortI[FldCntri].Sx; pY=VecSortI[FldCntri].Sy; pZ=VecSortI[FldCntri].Sz;
      if(VolLabl.value(pX, pY, pZ)==MASK){//sI
	curlab+=1;
	Qx.push(pX); Qy.push(pY); Qz.push(pZ);
	VolLabl.value(pX, pY, pZ)=curlab;
	while(!Qx.empty()){//sW
	  qX=Qx.front(); qY=Qy.front(); qZ=Qz.front();
	  Qx.pop(); Qy.pop(); Qz.pop();
	  counter=0;
	  for(int z0=-1; z0<=1; z0++)
	    for(int y0=-1; y0<=1; y0++)
	      for(int x0=-1; x0<=1; x0++){
		rX=qX+x0; rY=qY+y0; rZ=qZ+z0;
		doIT =(NumConn>=edsta[counter++]);                    
		if(doIT && (VolLabl.value(rX, rY, rZ)==MASK)){
		  Qx.push(rX); Qy.push(rY);Qz.push(rZ);
		  VolLabl.value(rX, rY, rZ)=curlab;
		}
	      }           
	}//eW       
      }//eI
    }
    ColumnVector ClusterSizes(curlab), ClusterSizesI(curlab);               
    ClusterSizes=0;
    for(tfceCntr=xFC; tfceCntr<sizeC; tfceCntr++){
      VecSortI[tfceCntr].Sl=VolLabl.value(VecSortI[tfceCntr].Sx, VecSortI[tfceCntr].Sy, VecSortI[tfceCntr].Sz);
      if ( VecSortI[tfceCntr].Sl>0 )
	ClusterSizes(VecSortI[tfceCntr].Sl)+=1;
    }
    float HH=pow(curThr, H);
    ClusterSizesI=pow(ClusterSizes, E)*HH;
    for(tfceCntr=xFC; tfceCntr<sizeC; tfceCntr++){
      if ( VecSortI[tfceCntr].Sl>0 ){
	VolEnhn.value(VecSortI[tfceCntr].Sx, VecSortI[tfceCntr].Sy, VecSortI[tfceCntr].Sz) += ClusterSizesI(VecSortI[tfceCntr].Sl);
      }
    }
  }//end curThr      
  copyconvert(VolEnhn,data);
  return;
}


template <class T>
void tfce_support(volume<T>& VolIntn, float H, float E, int NumConn, float minT, float deltaT, int Xoi, int Yoi, int Zoi, float threshTFCE)
{
  volume<float> VolEnhn; copyconvert(VolIntn, VolEnhn);
  volume<float> VolTemp; copyconvert(VolIntn, VolTemp);
  int minX=0, minY=0, minZ=0, maxX=VolIntn.maxx(), maxY=VolIntn.maxy(), maxZ=VolIntn.maxz();
  float maxT = VolIntn.value(Xoi,Yoi,Zoi);
  int thrCntr = 0; int thrNum = 0; 
  if(deltaT==0){
    deltaT = maxT/100;
    thrNum = 101;
  }
  else
    thrNum = int(maxT/deltaT)+1;
  ColumnVector Clusters(thrNum), Thresholds(thrNum), ClusterSizes;
  for(float thresh=minT; thresh<maxT; thresh+=deltaT){	
    copyconvert(VolEnhn, VolTemp);	
    VolTemp.binarise(thresh);
    volume<int> VolLabl = connected_components(VolTemp, ClusterSizes, NumConn);
    Clusters(thrCntr+1) = ClusterSizes(VolLabl.value(Xoi, Yoi, Zoi));
    Thresholds(thrCntr+1) = thresh;
    for(int z=minZ; z<=maxZ; z++)
      for(int y=minY; y<=maxY; y++)
	for(int x=minX; x<=maxX; x++)
	  if( VolLabl.value(x,y,z) != VolLabl.value(Xoi,Yoi,Zoi) )
	    VolEnhn.value(x,y,z) = 0;
    thrCntr++;
  }      
  float  summ=0, thre=0;
  for (int i=0; i<thrCntr; i++){
    summ+=pow(Clusters(thrCntr-i), E)*pow(Thresholds(thrCntr-i), H);
    thre = Thresholds(thrCntr-i);
    if(summ>=threshTFCE)
      break;
  }      
  if(summ<threshTFCE)
    cout<<"it doesn't reach to specified threshold"<<endl;
  copyconvert(VolIntn, VolEnhn); 
  copyconvert(VolIntn, VolTemp);     VolTemp.binarise(thre);
  volume<int> VolLabl = connected_components(VolTemp, ClusterSizes, NumConn);
  for(int z=minZ; z<=maxZ; z++)
    for(int y=minY; y<=maxY; y++)
      for(int x=minX; x<=maxX; x++)
	if( VolLabl.value(x,y,z)!=VolLabl(Xoi, Yoi, Zoi) )
	  VolEnhn.value(x,y,z) = 0;
  copyconvert(VolEnhn,VolIntn);
  return;
}
  



////////////////////////////////////////////////////////////////////////////
 	 
}
     

#endif


