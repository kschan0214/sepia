/*  newimage.cc

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 2000 University of Oxford  */

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

#include <complex>
#include <cassert>
#include <sstream>
#include <iostream>
#include "newmatio.h"
#include "newimage.h"
#include "fslio/fslio.h"
#include "arch.h"

using namespace NEWMAT;
using namespace LAZY;
using namespace MISCMATHS;

namespace NEWIMAGE {

  void imthrow(const string& msg, int nierrnum) {
    cerr << "Image Exception : #" << nierrnum << " :: " << msg << endl;
    if (USEASSERT) { bool no_error=false; assert(no_error); }
    else { throw Exception(msg.data()); }
  }

  // interfaces for lazy evaluated calculations
  template <class T>
  minmaxstuff<T> calc_minmax(const volume<T>& vol);

  template <class T>
  minmaxstuff<T> calc_minmax(const volume<T>& vol, const volume<T>& mask);

  template <class T>
  std::vector<double> calc_sums(const volume<T>& vol);

  template <class T>
  std::vector<double> calc_sums(const volume<T>& vol, const volume<T>& mask);

  template <class T>
  T calc_backgroundval(const volume<T>& vol);

  template <class T>
  ColumnVector calc_cog(const volume<T>& vol);

  template <class T>
  std::vector<T> calc_robustlimits(const volume<T>& vol);

  template <class T>
  Matrix calc_principleaxes(const volume<T>& vol);

  template <class T>
  std::vector<T> calc_percentiles(const volume<T>& vol);

  template <class T>
  std::vector<T> calc_percentiles(const volume<T>& vol, const volume<T>& mask, 
				  const std::vector<float>& percentilepvals);

  template <class T>
  ColumnVector calc_histogram(const volume<T>& vol);

  template <class T>
  ColumnVector calc_histogram(const volume<T>& vol, const volume<T>& mask);

  template <class T>
  ColumnVector calc_histogram(const volume4D<T>& vol);

  template <class T>
  ColumnVector calc_histogram(const volume4D<T>& vol, const volume4D<T>& mask);

  template <class T>
  ColumnVector calc_histogram(const volume4D<T>& vol, const volume<T>& mask);

  template<class T>
  SPLINTERPOLATOR::Splinterpolator<T> calc_spline_coefs(const volume<T>& vol);

  // Declaration of helper functions.

  SPLINTERPOLATOR::ExtrapolationType translate_extrapolation_type(extrapolation ep);
  
  // CONSTRUCTORS (not including copy constructor - see under copying)

  template <class T>
  int volume<T>::initialize(int xsize, int ysize, int zsize, T *d, bool d_owner)
    {
      // set up storage
      this->destroy();
      SlicesZ = zsize;
      RowsY = ysize;
      ColumnsX = xsize;
      SizeBound = SlicesZ * RowsY * ColumnsX;
      SliceOffset = RowsY * ColumnsX;
      if (SizeBound > 0) {
	if (d != 0) {
	  Data = d;
	  data_owner = d_owner;
	}
	else {
          try {
	    Data = new T[SizeBound];
	  } catch(...) { Data=0; }
	  if (Data==0) { imthrow("Out of memory",99); }
	  data_owner = true;
	}
      } else {
	Data = 0;
	data_owner = false;
      }
      
      setdefaultproperties();
      return 0;
    }


  template <class T>
  void volume<T>::setdefaultproperties()
    {
      Xdim = 1.0;
      Ydim = 1.0;
      Zdim = 1.0;

      StandardSpaceCoordMat = IdentityMatrix(4);
      RigidBodyCoordMat = IdentityMatrix(4);
      StandardSpaceTypeCode = NIFTI_XFORM_UNKNOWN;
      RigidBodyTypeCode = NIFTI_XFORM_UNKNOWN;
      RadiologicalFile = true;

      IntentCode = NIFTI_INTENT_NONE;
      IntentParam1 = 0.0;
      IntentParam2 = 0.0;
      IntentParam3 = 0.0;
      
      SliceOrderingCode = NIFTI_SLICE_UNKNOWN;
      
      Limits.resize(6,0);
      setdefaultlimits();
      // Default ROI is whole volume
      ROIbox = Limits;
      activeROI = false;
      calc_no_voxels();
      
      minmax.init(this,calc_minmax);
      sums.init(this,calc_sums);
      backgroundval.init(this,calc_backgroundval);
      lazycog.init(this,calc_cog);
      robustlimits.init(this,calc_robustlimits);
      principleaxes.init(this,calc_principleaxes);
      percentiles.init(this,calc_percentiles);
      l_histogram.init(this,calc_histogram);
      splint.init(this,calc_spline_coefs);
      HISTbins=256;
      HISTmin=(T) 0;
      HISTmax=(T) 0;

      // Initial percentile pvals to store when calculating percentiles
      percentilepvals.erase(percentilepvals.begin(),percentilepvals.end());
      percentilepvals.push_back(0.0);
      percentilepvals.push_back(0.001);
      percentilepvals.push_back(0.005);
      for (int probval=1; probval<=99; probval++) {
	percentilepvals.push_back(((float) probval)/100.0);
      }
      percentilepvals.push_back(0.995);
      percentilepvals.push_back(0.999);
      percentilepvals.push_back(1.0);

      p_interpmethod = trilinear;
      p_extrapmethod = zeropad;
      splineorder = 3;
      padvalue = (T) 0;
      extrapval = padvalue;
      p_userinterp = 0;
      p_userextrap = 0;
      ep_valid.resize(3);
      ep_valid[0] = false; ep_valid[1] = false; ep_valid[2] = false;

      displayMaximum=0;
      displayMinimum=0;
      strncpy(auxFile,string("").c_str(),24);

      set_whole_cache_validity(false);
    }


  template <class T>
  volume<T>::volume() : Data(0), data_owner(false)
    {
      this->initialize(0,0,0,0,false);
    }
  

  template <class T>
  volume<T>::volume(int xsize, int ysize, int zsize) 
    : Data(0), data_owner(false)
    {
      this->initialize(xsize,ysize,zsize,0,true);
    }


  template <class T>
  volume<T>::volume(int xsize, int ysize, int zsize, T *d, bool d_owner) 
    : Data(d), data_owner(false)
    {
      this->initialize(xsize,ysize,zsize,d,d_owner);
    }


  template <class T>
  int volume<T>::reinitialize(int xsize, int ysize, int zsize)
    {
      return this->initialize(xsize,ysize,zsize,0,true);
    }

  template <class T>
  int volume<T>::reinitialize(int xsize, int ysize, int zsize, T *d, 
			      bool d_owner)
    {
      return this->initialize(xsize,ysize,zsize,d,d_owner);
    }

  template <class T>
  void volume<T>::destroy()
  {
    if (data_owner) delete [] Data;
    Data = 0;
  }
  

  template <class T>
  volume<T>::~volume()
    {  this->destroy(); }



  template <class T>
  void volume<T>::calc_no_voxels() const 
  {
    no_voxels = ((unsigned long int)(maxx()-minx()+1)) *((unsigned long int) (maxy()-miny()+1)) *((unsigned long int) (maxz()-minz()+1)); 
  }
  
  template <class T>
  void volume<T>::setupsizeproperties() const 
  {
    calc_no_voxels(); 
  }

  template <class T>
  void volume<T>::setdefaultlimits() const
    {
      Limits[0]=0; Limits[1]=0; Limits[2]=0; 
      Limits[3]=this->xsize()-1; 
      Limits[4]=this->ysize()-1; 
      Limits[5]=this->zsize()-1;
    }

  
  template <class T>
  void volume<T>::enforcelimits(std::vector<int>& lims) const
    {
      // clamp all elements within valid bounds
//        lims[0]=Max(0,lims[0]); 
//        lims[1]=Max(0,lims[1]); 
//        lims[2]=Max(0,lims[2]); 
//        lims[3]=Min(this->xsize() - 1,lims[3]); 
//        lims[4]=Min(this->ysize() - 1,lims[4]); 
//        lims[5]=Min(this->zsize() - 1,lims[5]); 
    }

  // COPYING AND CONVERSION FUNCTIONS

  template <class T>
  int volume<T>::reinitialize(const volume<T>& source)
    {
      this->initialize(source.xsize(),source.ysize(),source.zsize(),0,false);
      this->copydata(source);
      return this->copyproperties(source);
    }
  
  template <class T>
  volume<T>::volume(const volume<T>& source) : Data(0), 
    data_owner(false)
    {
      this->reinitialize(source);
    }

  template <class T>
  const volume<T>& volume<T>::operator=(const volume<T>& source)
  {
    this->reinitialize(source);
    return *this;
  }

  template <class T>
  int volume<T>::copydata(const volume<T>& source) {
    if (SizeBound != source.SizeBound) {
      imthrow("Attempted to copydata with non-matching sizes",2);
      //return -1;
    }
    copy(source.Data, source.Data + SizeBound, Data);  // use the STL
    data_owner = true;
    return 0;
  }


  template <class T>
  volume<T> volume<T>::ROI() const
  {
    volume<T> roivol;
    roivol.reinitialize(this->maxx() - this->minx()+1,
			this->maxy() - this->miny()+1,
			this->maxz() - this->minz()+1,0,false);

    // now copy only the appropriate data
    for (int z=this->minz(); z<=this->maxz(); z++) {
      for (int y=this->miny(); y<=this->maxy(); y++) {
	for (int x=this->minx(); x<=this->maxx(); x++) {
	  roivol(x - this->minx(),y - this->miny(),z - this->minz()) =
	    (*this)(x,y,z);
	}
      }
    }
    roivol.copyproperties(*this);
    roivol.deactivateROI();
    // set sform and qform matrices appropriately (if set)
    Matrix roi2vol= IdentityMatrix(4);
    roi2vol(1,4) = this->minx();
    roi2vol(2,4) = this->miny();
    roi2vol(3,4) = this->minz();
    if (this->sform_code()!=NIFTI_XFORM_UNKNOWN) {
      roivol.set_sform(this->sform_code(),this->sform_mat() * roi2vol);
    }
    if (this->qform_code()!=NIFTI_XFORM_UNKNOWN) {
      roivol.set_qform(this->qform_code(),this->qform_mat() * roi2vol);
    }
    roivol.set_whole_cache_validity(false);
    return roivol;
  }


  template <class T>
  int volume<T>::copyROIonly(const volume<T>& source)
  {
    if (!samesize(*this,source)) {
      imthrow("Attempted to copy ROIs when different sizes",3);
    }

    // now copy only the appropriate data
    int xoff=source.minx()-minx(), yoff=source.miny()-miny(), 
      zoff=source.minz()-minz();
    for (int z=source.minz(); z<=source.maxz(); z++) {
      for (int y=source.miny(); y<=source.maxy(); y++) {
	for (int x=source.minx(); x<=source.maxx(); x++) {
	  (*this)(x - xoff,y - yoff,z - zoff) = source(x,y,z);
	}
      }
    }
    set_whole_cache_validity(false);
    return 0;
  }


  template <class T>
  int volume<T>::copyproperties(const volume<T>& source) 
  {
    // sets all properties
    copybasicproperties(source,*this);

    this->copylazymanager(source);
    minmax.copy(source.minmax,this);
    sums.copy(source.sums,this);
    backgroundval.copy(source.backgroundval,this);
    lazycog.copy(source.lazycog,this);
    robustlimits.copy(source.robustlimits,this);
    principleaxes.copy(source.principleaxes,this);
    percentiles.copy(source.percentiles,this);
    l_histogram.copy(source.l_histogram,this);
    HISTbins = source.HISTbins;
    HISTmin = source.HISTmin;
    HISTmax = source.HISTmax;    
    percentilepvals = source.percentilepvals;
    p_userextrap = source.p_userextrap;
    p_userinterp = source.p_userinterp;

    return 0;
  }

  // Volume->ColumnVector

  template<class T>
  ReturnMatrix volume<T>::vec(const volume<T>& mask) const
  {
    if (!samesize(mask,*this)) {imthrow("volume<T>::vec: Mask and volume of different size",3);}
    ColumnVector ovec(xsize()*ysize()*zsize());
    for (int vindx=0,k=0; k<zsize(); k++) {
      for (int j=0; j<ysize(); j++) {
	for (int i=0; i<xsize(); i++) {
          ovec.element(vindx) = (mask(i,j,k)>0) ? (*this)(i,j,k) : 0.0;
          vindx++;
	}
      }
    }
    ovec.Release();
    return ovec;
  }

  // Code multiplication to avoid allocating mask volume

  template <class T>
  ReturnMatrix volume<T>::vec() const
  {
    ColumnVector ovec(xsize()*ysize()*zsize());
    for (int vindx=0, k=0; k<zsize(); k++) {
      for (int j=0; j<ysize(); j++) {
	for (int i=0; i<xsize(); i++) {
          ovec.element(vindx) = (*this)(i,j,k);
          vindx++;
	}
      }
    }
    ovec.Release();
    return ovec;
  }

  // Insert ColumnVector into volume

  template <class T>
  void volume<T>::insert_vec(const ColumnVector&  pvec,
                             const volume<T>&     mask)
  {
    if (pvec.Nrows() != xsize()*ysize()*zsize()) {
      cout << "pvec.Nrows() = " << pvec.Nrows() << endl;
      cout << "xsize() = " << xsize() << ",  ysize() = " << ysize() << ",  zsize() = " << zsize() << endl;
      imthrow("volume<T>::insert_vec: Size mismatch between ColumnVector and image volume",3);
    }
    if (!samesize(mask,*this)) {
      imthrow("volume<T>::insert_vec: Size mismatch between mask and image volume",3);
    }
    for (int vindx=0, k=0; k<zsize(); k++) {
      for (int j=0; j<ysize(); j++) {
	for (int i=0; i<xsize(); i++) {
	  (*this)(i,j,k) = (mask(i,j,k) > 0) ? ((T) pvec.element(vindx)) : ((T) 0);
          vindx++;
	}
      }
    }
  }

  template <class T>
  void volume<T>::insert_vec(const ColumnVector&  pvec)
  {
    if (pvec.Nrows() != xsize()*ysize()*zsize()) {
      cout << "pvec.Nrows() = " << pvec.Nrows() << endl;
      cout << "xsize() = " << xsize() << ",  ysize() = " << ysize() << ",  zsize() = " << zsize() << endl;
      imthrow("volume<T>::insert_vec: Size mismatch between ColumnVector and image volume",3);
    }
    for (int vindx=0, k=0; k<zsize(); k++) {
      for (int j=0; j<ysize(); j++) {
	for (int i=0; i<xsize(); i++) {
	  (*this)(i,j,k) = ((T) pvec.element(vindx));
          vindx++;
	}
      }
    }
  }

  template <class T>
  vector<int> volume<T>::labelToCoord(const long label) const
  {
    vector<int> coordinates;
    coordinates.push_back(label%this->xsize());
    coordinates.push_back( (floor) ( ( label%( this->xsize()*this->ysize() ) ) / this->xsize() ));
    coordinates.push_back( (floor) ( label / ( this->xsize()*this->ysize() ) ) );
    return coordinates;
  }

  // ROI functions

  template <class T>
  void volume<T>::setROIlimits(int x0, int y0, int z0, 
				    int x1, int y1, int z1) const
    { 
      // Enforce ordering
      ROIbox[0]=Min(x0,x1); 
      ROIbox[1]=Min(y0,y1); 
      ROIbox[2]=Min(z0,z1); 
      ROIbox[3]=Max(x0,x1); 
      ROIbox[4]=Max(y0,y1); 
      ROIbox[5]=Max(z0,z1);
      enforcelimits(ROIbox);
      if (activeROI) activateROI();
    }

  template <class T>
  void volume<T>::setROIlimits(const std::vector<int>& lims) const
      { 
	if (lims.size()!=6) return; 
	setROIlimits(lims[0],lims[1],lims[2],lims[3],lims[4],lims[5]);
      }
 
  template <class T>
  void volume<T>::activateROI() const
  { 
    activeROI=true; 
    enforcelimits(ROIbox);
    Limits = ROIbox;
    set_whole_cache_validity(false);
    setupsizeproperties(); 
  }

  template <class T>
  void volume<T>::deactivateROI() const
  { 
    activeROI=false; 
    setdefaultlimits();
    set_whole_cache_validity(false); 
    setupsizeproperties(); 
  }


  int mirrorclamp(int x, int x1, int x2) {
    if (x2<x1) return mirrorclamp(x,x2,x1);
    if (x1==x2) return x1;
    int x3 = 2*x2 - x1 + 1;
    int nx = periodicclamp(x,x1,x3);
    if (nx > x2) {
      nx = 2*x2 + 1 - nx;
    }
    return nx;
  }

  // EXTRAPOLATION AND INTERPOLATION

  template <class T>
  void volume<T>::setinterpolationmethod(interpolation interpmethod) const
    { 
      p_interpmethod = interpmethod;
      // define a default sinc kernel if no kernel has previously been defined
      if ( (interpmethod == sinc) && (interpkernel.kernelvals()==0) ) {
	string sincwindowtype = "blackman";
	this->definesincinterpolation(sincwindowtype,7);
      }
    }

  template<class T>
  void volume<T>::setsplineorder(unsigned int order) const
  {
    if (order > 7) imthrow("setsplineorder: Only splines of order up to 7 allowed",10);
    splineorder = order;
  }

  template <class T>
  bool in_neigh_bounds(const volume<T>& vol, int x, int y, int z)
    {  return ( (x>=0) && (y>=0) && (z>=0) &&
		(x<(vol.xsize()-1)) && (y<(vol.ysize()-1)) && 
		(z<(vol.zsize()-1)) ); }

  inline float q_tri_interpolation(float v000, float v001, float v010, 
				   float v011, float v100, float v101, 
				   float v110, float v111, 
				   float dx, float dy, float dz)

    {
      float temp1, temp2, temp3, temp4, temp5, temp6;
      temp1 = (v100 - v000)*dx + v000;
      temp2 = (v101 - v001)*dx + v001;
      temp3 = (v110 - v010)*dx + v010;
      temp4 = (v111 - v011)*dx + v011;
      // second order terms
      temp5 = (temp3 - temp1)*dy + temp1;
      temp6 = (temp4 - temp2)*dy + temp2;
      // final third order term
      return (temp6 - temp5)*dz + temp5;
    }
  
  //////// Kernel Interpolation Call /////////
  
  template <class T>
  float volume<T>::kernelinterpolation(const float x, const float y, 
				       const float z) const
  {
    const kernelstorage* storedkernel = interpkernel.kernelvals();
    // sanity check on kernel
    if (storedkernel==0) {
      cerr << "ERROR: Must set kernel parameters before using interpolation!" 
	   << endl;
      return (float) extrapolate(0,0,0);
    }
    
    // kernel half-width  (i.e. range is +/- w)
    int wx=storedkernel->widthx();  
    int wy=storedkernel->widthy();  
    int wz=storedkernel->widthz();  
    ColumnVector kernelx = storedkernel->kernelx();
    ColumnVector kernely = storedkernel->kernely();
    ColumnVector kernelz = storedkernel->kernelz();
    float *storex = storedkernel->storex;
    float *storey = storedkernel->storey;
    float *storez = storedkernel->storez;
    
    int ix0, iy0, iz0;
    ix0 = (int) floor(x);
    iy0 = (int) floor(y);
    iz0 = (int) floor(z);
    
    float convsum=0.0, interpval=0.0, kersum=0.0;
    
    for (int d=-wz; d<=wz; d++) {
      storez[d+wz] = kernelval((z-iz0+d),wz,kernelz);
    }
    for (int d=-wy; d<=wy; d++) {
      storey[d+wy] = kernelval((y-iy0+d),wy,kernely);
    }
    for (int d=-wx; d<=wx; d++) {
      storex[d+wx] = kernelval((x-ix0+d),wx,kernelx);
    }
    
    int xj, yj, zj;
    for (int z1=iz0-wz; z1<=iz0+wz; z1++) {
      zj=iz0-z1+wz;
      for (int y1=iy0-wy; y1<=iy0+wy; y1++) {
	yj=iy0-y1+wy;
	for (int x1=ix0-wx; x1<=ix0+wx; x1++) {
	  if (in_bounds(x1,y1,z1)) {
	    xj=ix0-x1+wx;
	    float kerfac = storex[xj] * storey[yj] * storez[zj];
	    convsum += this->operator()(x1,y1,z1) * kerfac;
	    kersum += kerfac;
	  } 
	}
      }
    }
    
    if ( (fabs(kersum)>1e-9) ) {
      interpval = convsum / kersum;
    } else {
      interpval = (float) extrapolate(ix0,iy0,iz0);
    }
    return interpval;
    
  }
  
  // The following routines are used to obtain an interpolated intensity value and
  // either a selected partial derivative (dx, dy or dz) or all partial derivatives
  // at the same location. The routine returning all derivatives is useful for 
  // non-linear registration of one subject to another (or to an atlas) and the 
  // routines returning a single derivative are useful e.g. for distortion correction.
  // Puss J

  template <class T>
  float volume<T>::interp1partial(// Input
                                  float x, float y, float z,    // Co-ordinates to get value for
                                  int     dir,                  // Direction for partial, 0->x, 1->y, 2->z
                                  // Output
                                  float  *pderiv)               // Derivative returned here
  const
  {
    if (getinterpolationmethod() != trilinear && getinterpolationmethod() != spline) {
      imthrow("Derivatives only implemented for tri-linear and spline interpolation",10);
    }
    if (dir < 0 || dir > 2) {
      imthrow("Ivalid derivative direction",11);
    }
    if (getinterpolationmethod() == trilinear) {
      int ix = ((int) floor(x));
      int iy = ((int) floor(y));
      int iz = ((int) floor(z));
      float dx = x - ((float) ix);
      float dy = y - ((float) iy);
      float dz = z - ((float) iz);
      float v000, v001, v010, v011, v100, v101, v110, v111;
      if (!in_neigh_bounds(*this,ix,iy,iz)) {   // We'll have to do some extrapolation
	v000 = (float) this->operator()(ix,iy,iz);
	v001 = (float) this->operator()(ix,iy,iz+1);
	v010 = (float) this->operator()(ix,iy+1,iz);
	v011 = (float) this->operator()(ix,iy+1,iz+1);
	v100 = (float) this->operator()(ix+1,iy,iz);
	v101 = (float) this->operator()(ix+1,iy,iz+1);
	v110 = (float) this->operator()(ix+1,iy+1,iz);
	v111 = (float) this->operator()(ix+1,iy+1,iz+1);
      }
      else {
	T t000, t001, t010, t011, t100, t101, t110, t111;
	this->getneighbours(ix,iy,iz,t000,t001,t010,t011,t100,t101,t110,t111);
	v000 = ((float) t000); v001 = ((float) t001); v010 = ((float) t010);  
	v011 = ((float) t011); v100 = ((float) t100); v101 = ((float) t101);  
	v110 = ((float) t110); v111 = ((float) t111); 
      }
      // The (seemingly silly) code multiplication below is to
      // ensure that in no case does calculating one of the partials
      // neccessitate any calculation over and above just calculating
      // the interpolated value.
      float tmp11, tmp12, tmp13, tmp14;
      float tmp21, tmp22;
      if (dir == 0) {            // df/dx
	float onemdz = 1.0-dz;
	tmp11 = onemdz*v000 + dz*v001;
	tmp12 = onemdz*v010 + dz*v011;
	tmp13 = onemdz*v100 + dz*v101;
	tmp14 = onemdz*v110 + dz*v111;
	tmp21 = (1.0-dy)*tmp11 + dy*tmp12;
	tmp22 = (1.0-dy)*tmp13 + dy*tmp14;
	*pderiv = tmp22 - tmp21;
	return((1.0-dx)*tmp21 + dx*tmp22);
      }
      else if (dir == 1) {       // df/dy 
	float onemdz = 1.0-dz;
	tmp11 = onemdz*v000 + dz*v001;
	tmp12 = onemdz*v010 + dz*v011;
	tmp13 = onemdz*v100 + dz*v101;
	tmp14 = onemdz*v110 + dz*v111;
	tmp21 = (1.0-dx)*tmp11 + dx*tmp13;
	tmp22 = (1.0-dx)*tmp12 + dx*tmp14;
	*pderiv = tmp22 - tmp21;
	return((1.0-dy)*tmp21 + dy*tmp22);
      }
      else if (dir == 2) {       // df/dz
	float onemdy = 1.0-dy;
	tmp11 = onemdy*v000 + dy*v010;
	tmp12 = onemdy*v001 + dy*v011;
	tmp13 = onemdy*v100 + dy*v110;
	tmp14 = onemdy*v101 + dy*v111;
	tmp21 = (1.0-dx)*tmp11 + dx*tmp13;
	tmp22 = (1.0-dx)*tmp12 + dx*tmp14;
	*pderiv = tmp22 - tmp21;
	return((1.0-dz)*tmp21 + dz*tmp22);
      }
    }
    else if (getinterpolationmethod() == spline) {
      return(spline_interp1partial(x,y,z,dir,pderiv));
    }
    return(-1.0); // Should not be reached. Just to stop compiler from complaining.
  }

  template <class T>
  float volume<T>::interp3partial(// Input
                                  float x, float y, float z,              // Co-ordinates to get value for
                                   // Output
                                  float *dfdx, float *dfdy, float *dfdz)  // Partials
  const
  {
    if (getinterpolationmethod() != trilinear && getinterpolationmethod() != spline) {
      imthrow("interp3partial: Derivatives only implemented for tri-linear and spline interpolation",10);
    }
    if (getinterpolationmethod() == trilinear) {
      int ix = ((int) floor(x));
      int iy = ((int) floor(y));
      int iz = ((int) floor(z));
      float dx = x - ((float) ix);
      float dy = y - ((float) iy);
      float dz = z - ((float) iz);
      float v000, v001, v010, v011, v100, v101, v110, v111;
      if (!in_neigh_bounds(*this,ix,iy,iz)) {   // We'll have to do some extrapolation
	v000 = (float) this->operator()(ix,iy,iz);
	v001 = (float) this->operator()(ix,iy,iz+1);
	v010 = (float) this->operator()(ix,iy+1,iz);
	v011 = (float) this->operator()(ix,iy+1,iz+1);
	v100 = (float) this->operator()(ix+1,iy,iz);
	v101 = (float) this->operator()(ix+1,iy,iz+1);
	v110 = (float) this->operator()(ix+1,iy+1,iz);
	v111 = (float) this->operator()(ix+1,iy+1,iz+1);
      }
      else {
	T t000, t001, t010, t011, t100, t101, t110, t111;
	this->getneighbours(ix,iy,iz,t000,t001,t010,t011,t100,t101,t110,t111);
	v000 = ((float) t000); v001 = ((float) t001); v010 = ((float) t010);  
	v011 = ((float) t011); v100 = ((float) t100); v101 = ((float) t101);  
	v110 = ((float) t110); v111 = ((float) t111); 
      }
      //
      // And do linear interpolation with calculation of all partials
      //
      float onemdz = 1.0-dz;
      float onemdy = 1.0-dy;    
      float tmp11 = onemdz*v000 + dz*v001;
      float tmp12 = onemdz*v010 + dz*v011;
      float tmp13 = onemdz*v100 + dz*v101;
      float tmp14 = onemdz*v110 + dz*v111;
      *dfdx = onemdy*(tmp13-tmp11) + dy*(tmp14-tmp12);
      *dfdy = (1.0-dx)*(tmp12-tmp11) + dx*(tmp14-tmp13);
      tmp11 = onemdy*v000 + dy*v010;
      tmp12 = onemdy*v001 + dy*v011;
      tmp13 = onemdy*v100 + dy*v110;
      tmp14 = onemdy*v101 + dy*v111;
      float tmp21 = (1.0-dx)*tmp11 + dx*tmp13;
      float tmp22 = (1.0-dx)*tmp12 + dx*tmp14;
      *dfdz = tmp22 - tmp21;
      return(onemdz*tmp21 + dz*tmp22);
    }
    else if (getinterpolationmethod() == spline) {
      return(spline_interp3partial(x,y,z,dfdx,dfdy,dfdz));
    }
    return(0.0);  // To silence compiler.
  }

  template <class T>
  float volume<T>::spline_interp1partial(// Input
                                         float x, float y, float z,    // Co-ordinates to get value for
                                         int     dir,                  // Direction for partial, 0->x, 1->y, 2->z
                                         // Output
                                         float  *deriv)               // Derivative returned here
  const
  {
    if (!in_bounds(x,y,z)) {
      extrapolation ep = getextrapolationmethod();
      if (ep == boundsassert) { *deriv=0.0; assert(false); extrapval = padvalue; return(extrapval); }
      else if (ep == boundsexception) imthrow("splineinterpolate: Out of bounds",1);
      else if (ep == zeropad) { *deriv=0.0; extrapval = static_cast<T>(0.0); return(extrapval); }
      else if (ep == constpad) { *deriv=0.0; extrapval = padvalue; return(extrapval); } 
    }

    T         partial = static_cast<T>(0.0);
    float     rval = 0.0;
    const SPLINTERPOLATOR::Splinterpolator<T>&  sp = splint();
    if (getsplineorder() != sp.Order() || translate_extrapolation_type(getextrapolationmethod()) != sp.Extrapolation(0)) {
      const SPLINTERPOLATOR::Splinterpolator<T>& spp = splint.force_recalculation();
      rval = static_cast<float>(spp(x,y,z,dir,&partial));
    }
    else rval = static_cast<float>(sp(x,y,z,dir,&partial));
    *deriv = static_cast<float>(partial); 
    return(rval);
  }  

  template<class T>
  float volume<T>::spline_interp3partial(// Input
                                         float x, float y, float z,              // Co-ordinates to get value for
                                         // Output
                                         float *dfdx, float *dfdy, float *dfdz)  // Partials
  const
  {
    if (!in_bounds(x,y,z)) {
      extrapolation ep = getextrapolationmethod();
      if (ep == boundsassert) { *dfdx=0.0; *dfdy=0.0; *dfdz=0.0; assert(false); extrapval = padvalue; return(extrapval); }
      else if (ep == boundsexception) imthrow("splineinterpolate: Out of bounds",1);
      else if (ep == zeropad) { *dfdx=0.0; *dfdy=0.0; *dfdz=0.0; extrapval = static_cast<T>(0.0); return(extrapval); }
      else if (ep == constpad) { *dfdx=0.0; *dfdy=0.0; *dfdz=0.0; extrapval = padvalue; return(extrapval); } 
    }

    static std::vector<T>   partials(3,0);
    float                   rval = 0.0;
    const SPLINTERPOLATOR::Splinterpolator<T>&  sp = splint();
    if (getsplineorder() != sp.Order() || translate_extrapolation_type(getextrapolationmethod()) != sp.Extrapolation(0)) {
      const SPLINTERPOLATOR::Splinterpolator<T>& spp = splint.force_recalculation();
      rval = static_cast<float>(spp.ValAndDerivs(x,y,z,partials));
    }
    else rval = static_cast<float>(sp.ValAndDerivs(x,y,z,partials));
    *dfdx = static_cast<float>(partials[0]); 
    *dfdy = static_cast<float>(partials[1]); 
    *dfdz = static_cast<float>(partials[2]); 
    return(rval);
  }

  template<class T>
  float volume<T>::splineinterpolate(float x, float y, float z) const
  {
    extrapolation ep = getextrapolationmethod();

    if (!in_bounds(x,y,z)) {
      if (ep == boundsassert) { assert(false); extrapval = padvalue; return(extrapval); }
      else if (ep == boundsexception) imthrow("splineinterpolate: Out of bounds",1);
      else if (ep == zeropad) { extrapval = static_cast<T>(0.0); return(extrapval); }
      else if (ep == constpad) { extrapval = padvalue; return(extrapval); } 
    }
    if (ep == extraslice) if (!in_extraslice_bounds(x,y,z)) { extrapval = padvalue; return(extrapval); }

    const SPLINTERPOLATOR::Splinterpolator<T>&  sp = splint();
    if (getsplineorder() != sp.Order() || translate_extrapolation_type(ep) != sp.Extrapolation(0)) {
      const SPLINTERPOLATOR::Splinterpolator<T>& spp = splint.force_recalculation();
      return(static_cast<float>(spp(x,y,z)));
    }
    return(static_cast<float>(sp(x,y,z)));
  }

  template <class T>
  float volume<T>::interpolate(float x, float y, float z) const
    {
      int ix, iy, iz;
      switch (p_interpmethod) {
      case userinterpolation:
	if (p_userinterp == 0) {
	  imthrow("No user interpolation method set",7);
	} else {
	  return (*p_userinterp)(*this,x,y,z);
	}
      case nearestneighbour:
	ix=roundf(x); iy=roundf(y); iz=roundf(z);
	return this->operator()(ix,iy,iz);
      case trilinear:
	{
	  ix=(int) floor(x); iy=(int) floor(y); iz=(int) floor(z);
	  if (in_neigh_bounds(*this,ix,iy,iz)) return interpolatevalue(x,y,z);
	  float dx=x-ix, dy=y-iy, dz=z-iz;
	  float v000=0, v001=0, v010=0, v011=0, v100=0, v101=0, v110=0, v111=0;
	  v000 = (float) this->operator()(ix,iy,iz);
	  v001 = (float) this->operator()(ix,iy,iz+1);
	  v010 = (float) this->operator()(ix,iy+1,iz);
	  v011 = (float) this->operator()(ix,iy+1,iz+1);
	  v100 = (float) this->operator()(ix+1,iy,iz);
	  v101 = (float) this->operator()(ix+1,iy,iz+1);
	  v110 = (float) this->operator()(ix+1,iy+1,iz);
	  v111 = (float) this->operator()(ix+1,iy+1,iz+1);
	  return q_tri_interpolation(v000,v001,v010,v011,v100,v101,v110,v111,
				     dx,dy,dz);
	}
      case sinc:
      case userkernel:
	{
	  return kernelinterpolation(x,y,z);
	}
      case spline:
        {
          return(splineinterpolate(x,y,z));
	}
      default:
	imthrow("Invalid interpolation method",6);
      }
      return 0.0;  // Should never get to here
    }


  template <class T>
  float volume<T>::interpolatevalue(float x, float y, float z) const
    {
      int ix, iy, iz;
      switch (p_interpmethod) {
      case userinterpolation:
	if (p_userinterp == 0) {
	  imthrow("No user interpolation method set",7);
	} else {
	  return (*p_userinterp)(*this,x,y,z);
	}
      case nearestneighbour:
	ix=roundf(x); iy=roundf(y); iz=roundf(z);
	return value(ix,iy,iz);
      case trilinear:
	{
	  ix=(int) floor(x); iy=(int) floor(y); iz=(int) floor(z);
	  float dx=x-ix, dy=y-iy, dz=z-iz;
	  T t000=0, t001=0, t010=0, t011=0, t100=0, t101=0, t110=0, t111=0;
	  float v000, v001, v010, v011, v100, v101, v110, v111;
	  this->getneighbours(ix,iy,iz,t000,t001,t010,t011,t100,t101,t110,t111);
	  v000=(float) t000; v001=(float) t001; v010=(float) t010; 
	  v011=(float) t011; v100=(float) t100; v101=(float) t101; 
	  v110=(float) t110; v111=(float) t111; 
	  return q_tri_interpolation(v000,v001,v010,v011,v100,v101,v110,v111,
				     dx,dy,dz);
	}
      case sinc:
      case userkernel:
	{
	  return kernelinterpolation(x,y,z);
	}
      case spline:
        {
          return(splineinterpolate(x,y,z));
	}
      default:
	imthrow("Invalid interpolation method",6);
      }
      return 0.0;  // Should never get to here
    }


  template <class T>
  const T& volume<T>::extrapolate(int x, int y, int z) const
  {
    
    switch (getextrapolationmethod()) {
    case userextrapolation:
      if (p_userextrap == 0) {
        imthrow("No user extrapolation method set",7);
      } else {
        extrapval = (*p_userextrap)(*this,x,y,z);
        return extrapval;
      }
    case zeropad:
      extrapval = (T) 0;
      return extrapval;
    case constpad:
      extrapval = padvalue;
      return extrapval;
    default:
      ; // do nothing
    }
    int nx=x, ny=y, nz=z;
    switch (getextrapolationmethod()) {
    case periodic:
      nx = periodicclamp(x,Limits[0],Limits[3]);
      ny = periodicclamp(y,Limits[1],Limits[4]);
      nz = periodicclamp(z,Limits[2],Limits[5]);
      return value(nx,ny,nz);
    case mirror:
      nx = mirrorclamp(x,Limits[0],Limits[3]);
      ny = mirrorclamp(y,Limits[1],Limits[4]);
      nz = mirrorclamp(z,Limits[2],Limits[5]);
      return value(nx,ny,nz);
    case extraslice:
      if (nx==Limits[0]-1) { nx=Limits[0]; }
      else { if (nx==Limits[3]+1) nx=Limits[3]; }
      if (ny==Limits[1]-1) { ny=Limits[1]; }
      else { if (ny==Limits[4]+1) ny=Limits[4]; }
      if (nz==Limits[2]-1) { nz=Limits[2]; }
      else { if (nz==Limits[5]+1) nz=Limits[5]; }
      if (in_bounds(nx,ny,nz)) { return value(nx,ny,nz); }
      else { extrapval = padvalue; return extrapval; }
    case boundsexception:
      if (!in_bounds(x,y,z)) {
        ostringstream msg;
        msg << "Out of Bounds at ("<<x<<","<<y<<","<<z<<")";
        imthrow(msg.str(),1);
      } else {
        return extrapval; 
      }
    case boundsassert:
      assert(in_bounds(x,y,z));
      return extrapval;
    default:
      imthrow("Invalid extrapolation method",6);
    }

    return extrapval;
  }


  template <class T>
  void volume<T>::defineuserinterpolation(float (*interp)(
                         const volume<T>& , float, float, float)) const
  { 
    p_userinterp = interp;
  }


  template <class T>
  void volume<T>::defineuserextrapolation(T (*extrap)(
                         const volume<T>& , int, int, int)) const
  { 
    p_userextrap = extrap;
  }


  template <class T>
  void volume<T>::definekernelinterpolation(const ColumnVector& kx, 
					    const ColumnVector& ky,
					    const ColumnVector& kz, 
					    int wx, int wy, int wz) const
  {
    // takes full-widths and converts all to half-widths
    int hwx = (wx-1)/2;
    int hwy = (wy-1)/2;
    int hwz = (wz-1)/2;
    interpkernel.setkernel(kx,ky,kz,hwx,hwy,hwz);
  }

  template <class T>
  void volume<T>::definekernelinterpolation(const volume<T>& vol) const
  {
    // copying like this is safe
    interpkernel = vol.interpkernel;
  }

  // Support Functions
  
  template <class T>
  void volume<T>::definesincinterpolation(const string& sincwindowtype,
					  int w, int nstore) const
  {
    // full width
    this->definesincinterpolation(sincwindowtype,w,w,w,nstore);
  }

  template <class T>
  void volume<T>::definesincinterpolation(const string& sincwindowtype,
					  int wx, int wy, int wz, 
					  int nstore) const
  {
    // full widths
    if (nstore<1) nstore=1;
    ColumnVector kx, ky, kz;
    // calculate kernels
    kx = sinckernel1D(sincwindowtype,wx,nstore);
    ky = sinckernel1D(sincwindowtype,wy,nstore);
    kz = sinckernel1D(sincwindowtype,wz,nstore);
    
    this->definekernelinterpolation(kx,ky,kz,wx,wy,wz);
  }
  

  // PROPERTIES


  template <class T>
  Matrix volume<T>::sampling_mat() const
  {
    Matrix samp=IdentityMatrix(4);
    samp(1,1) = xdim();
    samp(2,2) = ydim();
    samp(3,3) = zdim();
    // NOTE: no origin information is contained in this matrix!
    return samp;
  }
  
  template <class T>  
  Matrix volume4D<T>::sampling_mat() const
  {
    return this->operator[](0).sampling_mat();
  }

  template <class T>
  void volume<T>::set_sform(int sform_code, const Matrix& snewmat) const
  {
    StandardSpaceTypeCode = sform_code;
    StandardSpaceCoordMat = snewmat;
  }


  template <class T>
  void volume<T>::set_qform(int qform_code, const Matrix& qnewmat) const
  {
    RigidBodyTypeCode = qform_code;
    RigidBodyCoordMat = qnewmat;
  }


  template <class T>
  Matrix volume4D<T>::sform_mat() const 
  {
    return this->operator[](0).sform_mat();
  }

  template <class T>
  int volume4D<T>::sform_code() const 
  {
    return this->operator[](0).sform_code();
  }


  template <class T>
  Matrix volume4D<T>::qform_mat() const 
  {
    return this->operator[](0).qform_mat();
  }

  template <class T>
  int volume4D<T>::qform_code() const 
  {
    return this->operator[](0).qform_code();
  }

  template <class T>
  void volume4D<T>::set_sform(int sform_code, const Matrix& snewmat) const 
  {
    for (int t=0; t<this->tsize(); t++) {
      vols[t].set_sform(sform_code,snewmat);
    }
  }

  template <class T>
  void volume4D<T>::set_qform(int qform_code, const Matrix& qnewmat) const 
  {
    for (int t=0; t<this->tsize(); t++) {
      vols[t].set_qform(qform_code,qnewmat);
    }
  }


  template <class T>
  float volume<T>::intent_param(int n) const
  {
    float retval=0;
    if (n==1) { retval = IntentParam1; }
    if (n==2) { retval = IntentParam2; }
    if (n==3) { retval = IntentParam3; }
    return retval;
  }


  template <class T>
  void volume<T>::set_intent(int intent_code, float p1, float p2, float p3) const
  {
    IntentCode = intent_code;
    IntentParam1 = p1;
    IntentParam2 = p2;
    IntentParam3 = p3;
  }


  template <class T>
  int volume4D<T>::intent_code() const 
  {
    return this->operator[](0).intent_code();
  }

  template <class T>
  float volume4D<T>::intent_param(int n) const 
  {
    return this->operator[](0).intent_param(n);
  }

  template <class T>
  void volume4D<T>::set_intent(int intent_code, float p1, float p2, float p3) 
    const 
  {
    for (int t=0; t<this->tsize(); t++) {
      vols[t].set_intent(intent_code,p1,p2,p3);
    }
  }



  template <class T>
  ColumnVector volume<T>::principleaxis(int n) const
  {
    Matrix tmp = principleaxes();
    ColumnVector res = tmp.SubMatrix(1,3,n,n);
    return res;
  }
  
  template <class T>
  Matrix volume<T>::principleaxes_mat() const
  {
    return principleaxes();
  }
  

  int pval_index_end() { return -1; }

  template <class T>
  int get_pval_index(const std::vector<T>& pvals, float p) 
  {
    int idx=0;
    while (idx < (int) pvals.size()) {
      // success if p is near pvals[idx] by a relative factor of 0.001 or less
      if ( fabs((p-pvals[idx])/Max(1e-5,Min(pvals[idx],1-pvals[idx]))) < 0.001 )
	return idx;
      else
	idx++;
    }
    return pval_index_end(); 
  }
    

  template <class T>
  T volume<T>::percentile(float pvalue, const volume<T>& mask) const
  {
    if ((pvalue>1.0) || (pvalue<0.0)) 
      { imthrow("Percentiles must be in the range [0.0,1.0]",4); }
    std::vector<float> pvaluevec;
    std::vector<T> retval;
    pvaluevec.push_back(pvalue);
    retval = calc_percentiles(*this,mask,pvaluevec);
    return retval[0];
  }


  template <class T>
  T volume<T>::percentile(float pvalue) const
  {
    if ((pvalue>1.0) || (pvalue<0.0)) 
      { imthrow("Percentiles must be in the range [0.0,1.0]",4); }
    int idx = get_pval_index(percentilepvals,pvalue);
    if (idx==pval_index_end()) {
      percentilepvals.push_back(pvalue);
      idx = percentilepvals.size() - 1;
      percentiles.force_recalculation();
    }
    assert((idx>=0) && (idx < (int) percentilepvals.size()));
    return percentiles()[idx];
  }


 
  template <class T>
  std::vector<T> percentile_vec(std::vector<T>& hist, 
				const std::vector<float>& percentilepvals)
  {
    unsigned int numbins = hist.size();
    if (numbins==0) {
      //cerr << "ERROR:: Empty image" << endl; // why call this an error? when the image is empty this should still work.....
      hist.push_back((T) 0);
      return hist;
    }
    
    sort(hist.begin(),hist.end());
    
    std::vector<T> outputvals(percentilepvals.size());
    for (unsigned int n=0; n<percentilepvals.size(); n++) {
      unsigned int percentile = 
	(unsigned int) (((float) numbins) * percentilepvals[n]);
      if (percentile<0)  percentile=0;
      if (percentile>=numbins)  percentile=numbins-1;
      outputvals[n] = hist[percentile];
    }
    return outputvals;
  }
  

  template <class T>
  std::vector<T> calc_percentiles(const volume<T>& vol, const volume<T>& mask, 
				  const std::vector<float>& percentilepvals)
  {
    if (!samesize(vol,mask)) {
      imthrow("mask and vol have different sizes in calc_percentiles",3);
    }
    std::vector<T> hist;
    for (int z=vol.minz(); z<=vol.maxz(); z++) {
      for (int y=vol.miny(); y<=vol.maxy(); y++) {
	for (int x=vol.minx(); x<=vol.maxx(); x++) {
	  if (mask(x,y,z)>0.5) hist.push_back(vol(x,y,z));
	}
      }
    }
    return percentile_vec(hist,percentilepvals);
  }      
  

  template <class T>
  std::vector<T> calc_percentiles(const volume<T>& vol)
  {
    unsigned int numbins = (unsigned int) vol.nvoxels();
    unsigned int hindx = 0;
    std::vector<T> hist(numbins);
    for (int z=vol.minz(); z<=vol.maxz(); z++) {
      for (int y=vol.miny(); y<=vol.maxy(); y++) {
	for (int x=vol.minx(); x<=vol.maxx(); x++) {
	  hist[hindx++] = vol(x,y,z);
	}
      }
    }
    return percentile_vec(hist,vol.percentilepvalues());
  }      
  

  template <class T>
  ColumnVector volume<T>::histogram(int nbins, T minval, T maxval) const
  {
    bool sameparams = true;
    if (HISTbins != nbins) {
      HISTbins = nbins;
      sameparams = false;
    }
    if (HISTmin != minval) {
      HISTmin = minval;
      sameparams = false;
    }
    if (HISTmax != maxval) {
      HISTmax = maxval;
      sameparams = false;
    }
    if (!sameparams) {
      l_histogram.force_recalculation();
    }
    return l_histogram();
  }


  template <class T>
  ColumnVector volume<T>::histogram(int nbins) const
  {
    return histogram(nbins,robustmin(),robustmax());
  }

  template <class T>
  ColumnVector volume<T>::histogram(int nbins, T minval, T maxval, 
				    const volume<T>& mask) const
  {
    ColumnVector hist;
    calc_histogram(*this,nbins,minval,maxval,hist,mask);
    return hist;
  }

  template <class T>
  ColumnVector volume<T>::histogram(int nbins, const volume<T>& mask) const
  {
    return histogram(nbins,robustmin(),robustmax(),mask);
  }


  template <class T>
  T volume<T>::min(const volume<T>& mask) const
  {
    minmaxstuff<T> retval;
    retval = calc_minmax(*this,mask);
    return retval.min;	
  }

  template <class T>
  T volume<T>::max(const volume<T>& mask) const
  {
    minmaxstuff<T> retval;
    retval = calc_minmax(*this,mask);
    return retval.max;	
  }

  template <class T>
  int volume<T>::mincoordx(const volume<T>& mask) const
  {
    minmaxstuff<T> retval;
    retval = calc_minmax(*this,mask);
    return retval.minx;	
  }

  template <class T>
  int volume<T>::mincoordy(const volume<T>& mask) const
  {
    minmaxstuff<T> retval;
    retval = calc_minmax(*this,mask);
    return retval.miny;	
  }

  template <class T>
  int volume<T>::mincoordz(const volume<T>& mask) const
  {
    minmaxstuff<T> retval;
    retval = calc_minmax(*this,mask);
    return retval.minz;	
  }

  template <class T>
  int volume<T>::maxcoordx(const volume<T>& mask) const
  {
    minmaxstuff<T> retval;
    retval = calc_minmax(*this,mask);
    return retval.maxx;	
  }

  template <class T>
  int volume<T>::maxcoordy(const volume<T>& mask) const
  {
    minmaxstuff<T> retval;
    retval = calc_minmax(*this,mask);
    return retval.maxy;	
  }

  template <class T>
  int volume<T>::maxcoordz(const volume<T>& mask) const
  {
    minmaxstuff<T> retval;
    retval = calc_minmax(*this,mask);
    return retval.maxz;	
  }

  template <class T>
  double volume<T>::mean(const volume<T>& mask) const
  { 
    return sum(mask)/(Max((double) no_mask_voxels(mask),1.0));
  }


  template <class T>
  double volume<T>::variance(const volume<T>& mask) const
  { 
    if (no_mask_voxels(mask)>0) {
      double n=(double) no_mask_voxels(mask);
      return (n/Max(1.0,n-1))*(sumsquares(mask)/n - mean(mask)*mean(mask));
    } else {
      cerr << "ERROR:: Empty mask image" << endl;
      return 0;
    }
  }


  template <class T>
  minmaxstuff<T> calc_minmax(const volume<T>& vol)
  {
    T newmin, newmax;
    int newminx=vol.minx(), newminy=vol.miny(), newminz=vol.minz(), 
        newmaxx=vol.minx(), newmaxy=vol.miny(), newmaxz=vol.minz();
    newmin = newmax = vol(vol.minx(),vol.miny(),vol.minz());
    for (int z=vol.minz(); z<=vol.maxz(); z++) {
      for (int y=vol.miny(); y<=vol.maxy(); y++) {
        for (int x=vol.minx(); x<=vol.maxx(); x++) {
          T val = vol(x,y,z);
          if (newmin>val) { newmin=val; newminx=x; newminy=y; newminz=z; }
          else if (val>newmax) { newmax=val;  newmaxx=x; newmaxy=y; newmaxz=z;}
        }
      }
    }
    minmaxstuff<T> newminmax;
    newminmax.min = newmin;
    newminmax.max = newmax;
    newminmax.minx = newminx;
    newminmax.miny = newminy;
    newminmax.minz = newminz;
    newminmax.mint = 0;
    newminmax.maxx = newmaxx;
    newminmax.maxy = newmaxy;
    newminmax.maxz = newmaxz;
    newminmax.maxt = 0;
    return newminmax;
  }



  template <class T>
  minmaxstuff<T> calc_minmax(const volume<T>& vol, const volume<T>& mask)
  {
    if (!samesize(vol,mask)) {
      imthrow("calc_minmax:: mask and volume must be the same size",4);
    }
    bool valid=false;
    T newmin, newmax;
    int newminx=vol.minx(), newminy=vol.miny(), newminz=vol.minz(), 
        newmaxx=vol.minx(), newmaxy=vol.miny(), newmaxz=vol.minz();
    newmin = newmax = vol(vol.minx(),vol.miny(),vol.minz());
    for (int z=vol.minz(); z<=vol.maxz(); z++) {
      for (int y=vol.miny(); y<=vol.maxy(); y++) {
        for (int x=vol.minx(); x<=vol.maxx(); x++) {
	  if (mask.value(x,y,z)>(T) 0.5) {
            T val = vol.value(x,y,z);
            if (!valid || (newmin>val)) { newmin=val; newminx=x; newminy=y; newminz=z; }
            if (!valid || (val>newmax)) { newmax=val;  newmaxx=x; newmaxy=y; newmaxz=z;}
	    valid=true;
	  }
        }
      }
    }
    minmaxstuff<T> newminmax;
    if (valid) {
      newminmax.min = newmin;
      newminmax.max = newmax;
      newminmax.minx = newminx;
      newminmax.miny = newminy;
      newminmax.minz = newminz;
      newminmax.mint = 0;
      newminmax.maxx = newmaxx;
      newminmax.maxy = newmaxy;
      newminmax.maxz = newmaxz;
      newminmax.maxt = 0;
    } else { 
	// invalid return type
      cerr << "ERROR:: Empty mask image" << endl;
      newminmax.min = 0;
      newminmax.max = 0;
      newminmax.minx = -1;
      newminmax.miny = -1;
      newminmax.minz = -1;
      newminmax.mint = -1;
      newminmax.maxx = -1;
      newminmax.maxy = -1;
      newminmax.maxz = -1;
      newminmax.maxt = -1;
    }
    return newminmax;
  }



  template <class T>
  std::vector<double> calc_sums(const volume<T>& vol)
  {
    double sum=0, sum2=0, totsum=0, totsum2=0;
    long int n=0, nlim;
    nlim = (long int) sqrt((double) vol.nvoxels());
    if (nlim<100000) nlim=100000;
    if (vol.usingROI()) {
      for (int z=vol.minz(); z<=vol.maxz(); z++) {
	for (int y=vol.miny(); y<=vol.maxy(); y++) {
	  for (int x=vol.minx(); x<=vol.maxx(); x++) {
	    T val = vol.value(x,y,z);
	    sum += val;
	    sum2 += val*val;
	    n++;
	    if (n>nlim) { n=0; totsum+=sum; totsum2+=sum2; sum=0; sum2=0; }
	  }
	}
      }
      totsum+=sum;
      totsum2+=sum2;
    } else {
      for (typename volume<T>::fast_const_iterator it=vol.fbegin(),
	       itend = vol.fend();   it!=itend; ++it) 
	{
	  T val = *it;
	  sum += val;
	  sum2 += val*val;
	  n++;
	  if (n>nlim) { n=0; totsum+=sum; totsum2+=sum2; sum=0; sum2=0; }
	}
      totsum+=sum;
      totsum2+=sum2;
    }
    std::vector<double> newsums(2);
    newsums[0] = totsum;
    newsums[1] = totsum2;
    return newsums;
  }


  template <class T>
  double volume<T>::sum(const volume<T>& mask) const
  {
    std::vector<double> retval;
    retval = calc_sums(*this,mask);
    return retval[0];	
  }

  template <class T>
  double volume<T>::sumsquares(const volume<T>& mask) const
  {
    std::vector<double> retval;
    retval = calc_sums(*this,mask);
    return retval[1];	
  }


  template <class T>
  std::vector<double> calc_sums(const volume<T>& vol, const volume<T>& mask)
  {
    if (!samesize(vol,mask)) {
      imthrow("calc_sums:: mask and volume must be the same size",4);
    }
    double sum=0, sum2=0, totsum=0, totsum2=0;
    long int n=0, nlim, nn=0;
    nlim = (long int) sqrt((double) vol.nvoxels());
    if (nlim<100000) nlim=100000;
    for (int z=vol.minz(); z<=vol.maxz(); z++) {
	for (int y=vol.miny(); y<=vol.maxy(); y++) {
	  for (int x=vol.minx(); x<=vol.maxx(); x++) {
	    if (mask.value(x,y,z)>(T) 0.5) {
	      T val = vol.value(x,y,z);
	      sum += val;
	      sum2 += val*val;
	      n++;
	      if (n>nlim) { nn++; n=0; totsum+=sum; totsum2+=sum2; sum=0; sum2=0; }
	    }
	  }
	}
      }
    totsum+=sum;
    totsum2+=sum2;
    std::vector<double> newsums(2);
    newsums[0] = totsum;
    newsums[1] = totsum2;
    if (n + nn == 0) {
      cerr << "ERROR:: Empty mask image" << endl;
    }
    return newsums;
  }


  template <class T>
  ColumnVector volume<T>::cog(const string& coordtype) const
  {
    // for coordtype="scaled_mm" return the old style, otherwise
    //  return newimage voxel coordinates
    ColumnVector retcog;
    retcog = this->lazycog();
    if (coordtype=="scaled_mm") {
      ColumnVector v(4);
      v << retcog(1) << retcog(2) << retcog(3) << 1.0;
      v = this->sampling_mat() * v;
      retcog(1) = v(1); retcog(2) = v(2); retcog(3) = v(3);
    }
    return retcog;
  }



  
  // the following calculates a robust background by taking the 10th percentile
  //  of the edge voxels only
  // Note: it does NOT use the ROI even if it is active
  template <class T>
  T calc_bval(const volume<T>& vol, unsigned int edgewidth)
    {
      unsigned int zb = vol.zsize(), yb = vol.ysize(), xb = vol.xsize();
      unsigned int ewx, ewy, ewz, numbins;
      ewx = edgewidth;  ewy = edgewidth;  ewz = edgewidth;
      if (ewx >= xb)  ewx=xb-1;
      if (ewy >= yb)  ewy=yb-1;
      if (ewz >= zb)  ewz=zb-1;
      numbins = 2*(xb-2*ewx)*(yb-2*ewy)*ewz + 2*(xb-2*ewx)*zb*ewy + 2*yb*zb*ewx;
      std::vector<T> hist(numbins);
      // put the edge voxel values into the histogram
      unsigned int hindx = 0;
      // put in the faces
      // xy faces (small lids of the box)
      for (unsigned int e=0; e<ewz; e++) {
	for (unsigned int x=ewx; x<xb-ewx; x++) {
	  for (unsigned int y=ewy; y<yb-ewy; y++) {
	    hist[hindx++] = vol.value(x,y,e);
	    hist[hindx++] = vol.value(x,y,zb-1-e);
	  }
	}
      }
      // xz faces (smallish edge faces)
      for (unsigned int e=0; e<ewy; e++) {
	for (unsigned int x=ewx; x<xb-ewx; x++) {
	  for (unsigned int z=0; z<zb; z++) {
	    hist[hindx++] = vol.value(x,e,z);
	    hist[hindx++] = vol.value(x,yb-1-e,z);
	  }
	}
      }
      // yz faces (large edge faces)
      for (unsigned int e=0; e<ewx; e++) {
	for (unsigned int y=0; y<yb; y++) {
	  for (unsigned int z=0; z<zb; z++) {
	    hist[hindx++] = vol.value(e,y,z);
	    hist[hindx++] = vol.value(xb-1-e,y,z);
	  }
	}
      }
      sort(hist.begin(),hist.end());
      unsigned int percentile10 = numbins / 10;
      T v10 = hist[percentile10];
      return v10;
    }

  template <class T>
  T calc_backgroundval(const volume<T>& vol)
  {
    return calc_bval(vol,2);
  }


  template <class T>
  ColumnVector calc_cog(const volume<T>& vol)
    {
      ColumnVector v_cog(3);
      v_cog(1)=0.0;
      v_cog(2)=0.0;
      v_cog(3)=0.0;
      double val=0, total=0, vx=0, vy=0, vz=0, tot=0;
      T vmin=vol.min();
      long int n=0, nlim;
      nlim = (long int) sqrt((double) vol.nvoxels());
      if (nlim<1000) nlim=1000;
      for (int z=vol.minz(); z<=vol.maxz(); z++) {
	for (int y=vol.miny(); y<=vol.maxy(); y++) {
	  for (int x=vol.minx(); x<=vol.maxx(); x++) {
	    val = (double) (vol(x,y,z) - vmin);
	    vx += val*x;
	    vy += val*y;
	    vz += val*z;
	    tot += val;
	    n++;
	    if (n>nlim) {
	      n=0; total+=tot; v_cog(1)+=vx; v_cog(2)+=vy; v_cog(3)+=vz; 
	      tot=0; vx=0; vy=0; vz=0;
	    }
	  }
	}
      }
      total+=tot; v_cog(1)+=vx; v_cog(2)+=vy; v_cog(3)+=vz;
      if (fabs(total) < 1e-5) {
	cerr << "WARNING::in calculating COG, total = 0.0" << endl;
	total = 1.0;
      }
      v_cog(1) /= total;
      v_cog(2) /= total;
      v_cog(3) /= total;
      // Leave these values in (newimage) voxel coordinates
      return v_cog;
    }


  template <class T>
  Matrix calc_principleaxes(const volume<T>& vol)
    {
      SymmetricMatrix m2(3);
      m2 = 0;
      double val=0, total=0, tot=0;
      double mxx=0, mxy=0, mxz=0, myy=0, myz=0, mzz=0, mx=0, my=0, mz=0;
      ColumnVector mean(3);
      mean = 0;
      T vmin=vol.min();

      long int n=0, nlim;
      nlim = (long int) sqrt((double) vol.nvoxels());
      if (nlim<1000) nlim=1000;
      for (int z=vol.minz(); z<=vol.maxz(); z++) {
	for (int y=vol.miny(); y<=vol.maxy(); y++) {
	  for (int x=vol.minx(); x<=vol.maxx(); x++) {
	    val = (double) (vol(x,y,z) - vmin);
	    mxx += val*x*x;
	    mxy += val*x*y;
	    mxz += val*x*z;
	    myy += val*y*y;
	    myz += val*y*z;
	    mzz += val*z*z;
	    mx += val*x;
	    my += val*y;
	    mz += val*z;
	    tot += val;
	    n++;
	    if (n>nlim) {
	      n=0; total+=tot; m2(1,1)+=mxx; m2(1,2)+=mxy; m2(1,3)+=mxz;
	      m2(2,2)+=myy; m2(2,3)+=myz; m2(3,3)+=mzz;
	      mean(1)+=mx; mean(2)+=my; mean(3)+=mz;
	      tot=0; mxx=0; mxy=0; mxz=0; myy=0; myz=0; mzz=0; mx=0; my=0; mz=0;
	    }
	  }
	}
      }
      total+=tot; m2(1,1)+=mxx; m2(1,2)+=mxy; m2(1,3)+=mxz;
      m2(2,2)+=myy; m2(2,3)+=myz; m2(3,3)+=mzz;
      mean(1)+=mx; mean(2)+=my; mean(3)+=mz;

      if (fabs(total) < 1e-5) {
	cerr << "WARNING::in calculating Principle Axes, total = 0.0" << endl;
	total = 1.0;
      }
      m2 /= total;
      mean /= total;
      // Now adjust for voxel dimensions
      Matrix samp(3,3);
      samp = vol.sampling_mat().SubMatrix(1,3,1,3);
      m2 << samp * m2 * samp;
      mean = samp*mean;
      // Now make it central (taking off the cog)
      Matrix meanprod(3,3);
      for (int k1=1; k1<=3; k1++) {
	for (int k2=1; k2<=3; k2++) {
	  meanprod(k1,k2) = mean(k1)*mean(k2);
	}
      }
      m2 << m2 - meanprod;

      Matrix paxes;
      DiagonalMatrix evals;
      Jacobi(m2,evals,paxes);
      // Force the eigenvalues (and vectors) to be in descending order
      ColumnVector ptemp;
      float etemp;
      // brute force sort the eigen values and vectors
      // find inedx of least e-value
      int kmin=1;
      for (int k=2; k<=3; k++) {
	if (evals(k,k) < evals(kmin,kmin)) kmin = k;
      }
      // put the least in position 1
      etemp = evals(1,1);
      ptemp = paxes.SubMatrix(1,3,1,1);
      evals(1,1) = evals(kmin,kmin);
      paxes.SubMatrix(1,3,1,1) = paxes.SubMatrix(1,3,kmin,kmin);
      evals(kmin,kmin) = etemp;
      paxes.SubMatrix(1,3,kmin,kmin) = ptemp;
      // check if remaining ones require swapping
      if (evals(3,3) < evals(2,2)) {
	etemp = evals(2,2);
	ptemp = paxes.SubMatrix(1,3,2,2);
	evals(2,2) = evals(3,3);
	paxes.SubMatrix(1,3,2,2) = paxes.SubMatrix(1,3,3,3);
	evals(3,3) = etemp;
	paxes.SubMatrix(1,3,3,3) = ptemp;
      }
      return paxes;
    }




  template <class T>
  int calc_histogram(const volume<T>& vol, int nbins, double minval, 
		     double maxval, ColumnVector& hist, const volume<T>& mask, 
		     bool use_mask=true)
  {
    // MJ NOTE: Concerned about the behaviour for integer volumes
    //           as rounding could cause values to go into the wrong bins
    if (hist.Nrows()!=nbins) hist.ReSize(nbins);
    hist = 0.0;
    
    if (maxval < minval) return -1;
    
    double a=((double) nbins)/(maxval - minval);
    double b= - ((double) nbins)*minval/(maxval - minval);
    int binno = 0;
    
    for (int z=vol.minz(); z<=vol.maxz(); z++) { 
      for (int y=vol.miny(); y<=vol.maxy(); y++) { 
	for (int x=vol.minx(); x<=vol.maxx(); x++) {
	  if ( (!use_mask) || (mask(x,y,z)>(T) 0.5) ) {
	    binno = (int) (a*((double) vol(x,y,z)) + b);
	    if (binno > nbins-1)  binno=nbins-1;
	    if (binno < 0)      binno=0;
	    hist(binno+1)++;
	  }
	}
      }
    }
    return 0; 
  }
  



  template <class T>
  int calc_histogram(const volume<T>& vol, int nbins, double minval, 
		     double maxval, ColumnVector& hist)
  {
    return calc_histogram(vol,nbins,minval,maxval,hist,vol,false);
  }

  

  template <class T>
  int calc_histogram(const volume4D<T>& vol, int nbins, double minval, 
		     double maxval, ColumnVector& hist, const volume4D<T>& mask,
		     bool use_mask=true)
  {
    // MJ NOTE: Concerned about the behaviour for integer volumes
    //           as rounding could cause values to go into the wrong bins
    if (!samesize(vol[0],mask[0])) {
      imthrow("calc_histogram:: mask and volume must be the same size",4);
    }

    if (hist.Nrows()!=nbins) hist.ReSize(nbins);
    hist = 0.0;
    
    if (maxval < minval) return -1;
    
    double a=((double) nbins)/(maxval - minval);
    double b= - ((double) nbins)*minval/(maxval - minval);
    int binno = 0;
    
    for (int t=vol.mint(); t<=vol.maxt(); t++) { 
      for (int z=vol.minz(); z<=vol.maxz(); z++) { 
	for (int y=vol.miny(); y<=vol.maxy(); y++) { 
	  for (int x=vol.minx(); x<=vol.maxx(); x++) {
	    if ( (!use_mask) || (mask(x,y,z,Min(t,mask.maxt()))>(T) 0.5) ) {
	      binno = (int) (a*((double) vol(x,y,z,t)) + b);
	      if (binno > nbins-1)  binno=nbins-1;
	      if (binno < 0)      binno=0;
	      hist(binno+1)++;
	    }
	  }
	}
      }
    }
    return 0; 
  }
  


  template <class T>
  int calc_histogram(const volume4D<T>& vol, int nbins, double minval, 
		     double maxval, ColumnVector& hist, const volume<T>& mask,
		     bool use_mask=true)
  {
    // MJ NOTE: Concerned about the behaviour for integer volumes
    //           as rounding could cause values to go into the wrong bins
    if (!samesize(vol[0],mask)) {
      imthrow("calc_histogram:: mask and volume must be the same size",4);
    }

    if (hist.Nrows()!=nbins) hist.ReSize(nbins);
    hist = 0.0;
    
    if (maxval < minval) return -1;
    
    double a=((double) nbins)/(maxval - minval);
    double b= - ((double) nbins)*minval/(maxval - minval);
    int binno = 0;
    
    for (int t=vol.mint(); t<=vol.maxt(); t++) { 
      for (int z=vol.minz(); z<=vol.maxz(); z++) { 
	for (int y=vol.miny(); y<=vol.maxy(); y++) { 
	  for (int x=vol.minx(); x<=vol.maxx(); x++) {
	    if ( (!use_mask) || (mask(x,y,z)>(T) 0.5) ) {
	      binno = (int) (a*((double) vol(x,y,z,t)) + b);
	      if (binno > nbins-1)  binno=nbins-1;
	      if (binno < 0)      binno=0;
	      hist(binno+1)++;
	    }
	  }
	}
      }
    }
    return 0; 
  }
  

  template <class T>
  int calc_histogram(const volume4D<T>& vol, int nbins, double minval, 
		     double maxval, ColumnVector& hist)
  {
    return calc_histogram(vol,nbins,minval,maxval,hist,vol,false);
  }



  template <class T>
  int find_histogram(const volume<T>& vol, ColumnVector& hist, int bins, 
		     T& min, T& max)
  {
    // STEVE SMITH'S CODE - ADAPTED FOR NEWIMAGE BY MARK JENKINSON
    int validsize=0;
    /* zero histogram */
    hist=0;
    if (min==max) return -1;
    
    /* create histogram; the MIN is so that the maximum value falls in the last valid bin, not the (last+1) bin */

    double fA = ((double)bins)/(max-min);
    double fB = ( ((double)bins) * ((double)(-min)) ) / (max-min);
    for (int z=vol.minz(); z<=vol.maxz(); z++) {
      for (int y=vol.miny(); y<=vol.maxy(); y++) {
	for (int x=vol.minx(); x<=vol.maxx(); x++) {
	  hist(Max(0, Min( (int)(fA*(vol(x,y,z)) + fB), bins-1) ) + 1)++;
	  validsize++;
	}
      }
    }
    
    return validsize;
  }

  
  
  template <class T>
  int find_histogram(const volume<T>& vol, ColumnVector& hist, int bins, 
		     T& min, T& max, const volume<T>& mask)
  {
    // STEVE SMITH'S CODE - ADAPTED FOR NEWIMAGE BY MARK JENKINSON
    if (!samesize(vol,mask)) { 
      imthrow("find_histogram:: mask and volume must be the same size",4);
    }
    if (no_mask_voxels(mask)==0) {
      cerr << "ERROR:: Empty mask image" << endl;
      return 0;
    }
    int validsize=0;
    /* zero histogram */
    hist=0;
    if (min==max) return -1;
    
    /* create histogram; the MIN is so that the maximum value falls in the last valid bin, not the (last+1) bin */

    double fA = ((double)bins)/(max-min);
    double fB = ( ((double)bins) * ((double)(-min)) ) / (max-min);
    for (int z=vol.minz(); z<=vol.maxz(); z++) {
      for (int y=vol.miny(); y<=vol.maxy(); y++) {
	for (int x=vol.minx(); x<=vol.maxx(); x++) {
	  if ( (mask(x,y,z)>(T) 0.5) )
	    {
	      hist(Max(0, Min( (int)(fA*(vol(x,y,z)) + fB), bins-1) ) + 1)++;
	      validsize++;
	    }
	}
      }
    }
    
    return validsize;
  }
  



  template <class T>
  int find_histogram(const volume4D<T>& vol, ColumnVector& hist, int bins, 
		     T& min, T& max, const volume4D<T>& mask)
  {
    // STEVE SMITH'S CODE - ADAPTED FOR NEWIMAGE BY MARK JENKINSON
    if (!samesize(vol[0],mask[0])) { 
      imthrow("find_histogram:: mask and volume must be the same size",4);
    }
    if (no_mask_voxels(mask)==0) {
      cerr << "ERROR:: Empty mask image" << endl;
      return 0;
    }
    int validsize=0;
    /* zero histogram */
    hist=0;
    if (min==max) return -1;
    
    /* create histogram; the MIN is so that the maximum value falls in the last valid bin, not the (last+1) bin */

    double fA = ((double)bins)/(max-min);
    double fB = ( ((double)bins) * ((double)(-min)) ) / (max-min);
    for (int t=vol.mint(); t<=vol.maxt(); t++) {
      for (int z=vol.minz(); z<=vol.maxz(); z++) {
	for (int y=vol.miny(); y<=vol.maxy(); y++) {
	  for (int x=vol.minx(); x<=vol.maxx(); x++) {
	    if ( (mask(x,y,z,Min(t,mask.maxt()))>(T) 0.5) )
	      {
		hist(Max(0, Min( (int)(fA*(vol(x,y,z,t)) + fB), bins-1) ) + 1)++;
		validsize++;
	      }
	  }
	}
      }
    }
    
    return validsize;
  }
  


  template <class T>
  int find_histogram(const volume4D<T>& vol, ColumnVector& hist, int bins, 
		     T& min, T& max, const volume<T>& mask)
  {
    // STEVE SMITH'S CODE - ADAPTED FOR NEWIMAGE BY MARK JENKINSON
    if (!samesize(vol[0],mask)) { 
      imthrow("find_histogram:: mask and volume must be the same size",4);
    }
    if (no_mask_voxels(mask)==0) {
      cerr << "ERROR:: Empty mask image" << endl;
      return 0;
    }
    int validsize=0;
    /* zero histogram */
    hist=0;
    if (min==max) return -1;
    /* create histogram; the MIN is so that the maximum value falls in the last valid bin, not the (last+1) bin */
    double fA = ((double)bins)/(max-min);
    double fB = ( ((double)bins) * ((double)(-min)) ) / (max-min);
    for (int t=vol.mint(); t<=vol.maxt(); t++) {
      for (int z=vol.minz(); z<=vol.maxz(); z++) {
	for (int y=vol.miny(); y<=vol.maxy(); y++) {
	  for (int x=vol.minx(); x<=vol.maxx(); x++) {
	    if ( (mask(x,y,z)>(T) 0.5) ) 
	      {
		hist(Max(0, Min( (int)(fA*(vol(x,y,z,t)) + fB), bins-1) ) + 1)++;
		validsize++;
	      }
	  }
	}
      }
    }
    
    return validsize;
  }
  

  template <class T>
  int find_histogram(const volume4D<T>& vol, ColumnVector& hist, int bins, 
		     T& min, T& max)
  {
    // STEVE SMITH'S CODE - ADAPTED FOR NEWIMAGE BY MARK JENKINSON
    int validsize=0;
    /* zero histogram */
    hist=0;
    if (min==max) return -1;
    /* create histogram; the MIN is so that the maximum value falls in the last valid bin, not the (last+1) bin */

    double fA = ((double)bins)/(max-min);
    double fB = ( ((double)bins) * ((double)(-min)) ) / (max-min);
    for (int t=vol.mint(); t<=vol.maxt(); t++) {
      for (int z=vol.minz(); z<=vol.maxz(); z++) {
	for (int y=vol.miny(); y<=vol.maxy(); y++) {
	  for (int x=vol.minx(); x<=vol.maxx(); x++) {
	    hist(Max(0, Min( (int)(fA*(vol(x,y,z,t)) + fB), bins-1) ) + 1)++;
	    validsize++;
	  }
	}
      }
    }
    
    return validsize;
  }
  



  template <class T, class S, class R>
  void find_thresholds(const S& vol, T& minval, T& maxval, const R& mask, bool use_mask=true)
  {
    // STEVE SMITH'S CODE - ADAPTED FOR NEWIMAGE BY MARK JENKINSON
  int HISTOGRAM_BINS=1000; 
  ColumnVector hist(HISTOGRAM_BINS);
  int MAX_PASSES=10;
  int top_bin=0, bottom_bin=0, count, pass=1,
    lowest_bin=0, highest_bin=HISTOGRAM_BINS-1, validsize;
  T thresh98=0, thresh2=0, min, max;
  if (use_mask) { min=vol.min(mask), max=vol.max(mask); }
  else { min=vol.min();  max=vol.max(); }
  if (hist.Nrows()!=HISTOGRAM_BINS) { hist.ReSize(HISTOGRAM_BINS); }

  while ( (pass==1) ||
	  ( (double) (thresh98 - thresh2) < (((double) (max - min)) / 10.0) ) ) // test for very long tails
    // find histogram and thresholds
    { 
      if (pass>1) // redo histogram with new min and max
	{
	  // increase range slightly from the 2-98% range found
	  bottom_bin=Max(bottom_bin-1,0);           
	  top_bin=Min(top_bin+1,HISTOGRAM_BINS-1);
	  
	  // now set new min and max on the basis of this new range
	  T tmpmin = (T)( min + ((double)bottom_bin/(double)(HISTOGRAM_BINS))*(max-min) );
	  max = (T)( min + ((double)(top_bin+1)/(double)(HISTOGRAM_BINS))*(max-min) );
	  min=tmpmin;
	}

      if (pass==MAX_PASSES || min==max)  // give up and revert to full range ...
	{
	  if (use_mask) { min=vol.min(mask); max=vol.max(mask); }
 	  else { min=vol.min();  max=vol.max(); }
	}
       
      if (use_mask) validsize = find_histogram(vol,hist,HISTOGRAM_BINS,min,max,mask);
      else validsize = find_histogram(vol,hist,HISTOGRAM_BINS,min,max);
      
      if (validsize<1)
	{
          minval=thresh2=min;
	  maxval=thresh98=max;
	  return;
	}    
      
      if (pass==MAX_PASSES)  /* ... _but_ ignore end bins */
	{
	  validsize-= MISCMATHS::round(hist(lowest_bin+1)) + 
	    MISCMATHS::round(hist(highest_bin+1));
	  lowest_bin++;
	  highest_bin--;
	}
      
      if (validsize<0) /* ie zero range */
	{

	  thresh2=thresh98=min;
	  break;
	}
      
      double fA = (max-min)/(double)(HISTOGRAM_BINS);

      for(count=0, bottom_bin=lowest_bin; count<validsize/50; bottom_bin++)
	count+=MISCMATHS::round(hist(bottom_bin + 1));
      bottom_bin--;
      thresh2 =  min + (T)((double)bottom_bin*fA);

      for(count=0, top_bin=highest_bin; count<validsize/50; top_bin--)
	count+=MISCMATHS::round(hist(top_bin + 1));
      top_bin++;
      thresh98 = min + (T)((double)(top_bin+1)*fA);

      if (pass==MAX_PASSES) break;
      pass++;
    }
    minval=thresh2;
    maxval=thresh98;
   
  }


  template <class T, class S>
  void find_thresholds(const S& vol, T& minval, T& maxval) {
    return find_thresholds(vol,minval,maxval,vol,false);
  }



  template <class T>
  std::vector<T> calc_robustlimits(const volume<T>& vol)
  {
    std::vector<T> rlimits(2);
    T minval=0, maxval=0;
    find_thresholds(vol,minval,maxval);
    //find_robust_limits(vol,1000,hist,minval,maxval);  // MJ version
    rlimits[0] = minval;
    rlimits[1] = maxval;
    return rlimits;
  }

  template <class T>
  std::vector<T> calc_robustlimits(const volume<T>& vol, const volume<T>& mask)
  {
    std::vector<T> rlimits(2);
    if (no_mask_voxels(mask)==0) {
      cerr << "ERROR:: Empty mask image" << endl;
      rlimits[0]=0;
      rlimits[1]=0;
      return rlimits;
    }
    T minval=0, maxval=0;
    find_thresholds(vol,minval,maxval,mask);
    rlimits[0] = minval;
    rlimits[1] = maxval;
    return rlimits;
  }


  template <class T>
  std::vector<T> calc_robustlimits(const volume4D<T>& vol)
  {
    std::vector<T> rlimits(2);
    T minval=0, maxval=0;
    find_thresholds(vol,minval,maxval);
    rlimits[0] = minval;
    rlimits[1] = maxval;
    return rlimits;
  }


  template <class T>
  std::vector<T> calc_robustlimits(const volume4D<T>& vol, const volume4D<T>& mask)
  {
    std::vector<T> rlimits(2);
    if (no_mask_voxels(mask)==0) {
      cerr << "ERROR:: Empty mask image" << endl;
      rlimits[0]=0;
      rlimits[1]=0;
      return rlimits;
    }
    T minval=0, maxval=0;
    find_thresholds(vol,minval,maxval,mask);
    rlimits[0] = minval;
    rlimits[1] = maxval;
    return rlimits;
  }

  

  template <class T>
  std::vector<T> calc_robustlimits(const volume4D<T>& vol, const volume<T>& mask)
  {
    std::vector<T> rlimits(2);
    if (no_mask_voxels(mask)==0) {
      cerr << "ERROR:: Empty mask image" << endl;
      rlimits[0]=0;
      rlimits[1]=0;
      return rlimits;
    }
    T minval=0, maxval=0;
    find_thresholds(vol,minval,maxval,mask);
    rlimits[0] = minval;
    rlimits[1] = maxval;
    return rlimits;
  }

  

  template <class T>
  T volume<T>::robustmin(const volume<T>& mask) const
  {
    std::vector<T> rlim;
    rlim = calc_robustlimits(*this,mask);
    return rlim[0];
  }

  template <class T>
  T volume<T>::robustmax(const volume<T>& mask) const
  {
    std::vector<T> rlim;
    rlim = calc_robustlimits(*this,mask);
    return rlim[1];
  }

  template <class T>
  T volume4D<T>::robustmin(const volume<T>& mask) const
  {
    std::vector<T> rlim;
    rlim = calc_robustlimits(*this,mask);
    return rlim[0];
  }

  template <class T>
  T volume4D<T>::robustmax(const volume<T>& mask) const
  {
    std::vector<T> rlim;
    rlim = calc_robustlimits(*this,mask);
    return rlim[1];
  }


  template <class T>
  T volume4D<T>::robustmin(const volume4D<T>& mask) const
  {
    std::vector<T> rlim;
    rlim = calc_robustlimits(*this,mask);
    return rlim[0];
  }

  template <class T>
  T volume4D<T>::robustmax(const volume4D<T>& mask) const
  {
    std::vector<T> rlim;
    rlim = calc_robustlimits(*this,mask);
    return rlim[1];
  }


  template <class T>
  ColumnVector calc_histogram(const volume<T>& vol)
  {
    ColumnVector hist;
    calc_histogram(vol,vol.histbins(),vol.histmin(),vol.histmax(),hist);
    return hist;
  }


  template <class T>
  ColumnVector calc_histogram(const volume<T>& vol, const volume<T>& mask)
  {
    ColumnVector hist;
    calc_histogram(vol,vol.histbins(),vol.histmin(),vol.histmax(),hist,mask);
    return hist;
  }

  template<class T>
  SPLINTERPOLATOR::Splinterpolator<T> calc_spline_coefs(const volume<T>& vol)
  {
    std::vector<unsigned int>                        dim(3,0);
    dim[0] = vol.xsize(); dim[1] = vol.ysize(); dim[2] = vol.zsize();
    std::vector<SPLINTERPOLATOR::ExtrapolationType>  ep(3,SPLINTERPOLATOR::Mirror);
    for (unsigned int i=0; i<3; i++) ep[i] = translate_extrapolation_type(vol.getextrapolationmethod());
    
    SPLINTERPOLATOR::Splinterpolator<T>  rval(vol.fbegin(),dim,ep,vol.getsplineorder(),false);

    return(rval);
  }

  SPLINTERPOLATOR::ExtrapolationType translate_extrapolation_type(extrapolation ep)
  {
    switch (ep) {
    case zeropad:
      return(SPLINTERPOLATOR::Zeros);
      break;
    case extraslice:
      return(SPLINTERPOLATOR::Constant);  // It is constant for a given column, hence name.
      break;
    case mirror:
      return(SPLINTERPOLATOR::Mirror);
      break;
    case periodic:
      return(SPLINTERPOLATOR::Periodic);
      break;
    case boundsassert: case boundsexception: // We deal with this at the actual interpolation, and for now just return something
      return(SPLINTERPOLATOR::Zeros);
      break;
    case constpad: // Not implemented in splinterpolator, so I'll deal with this too at the actual interpolation.
      return(SPLINTERPOLATOR::Zeros); 
      break;
    case userextrapolation:
      imthrow("translate_extrapolation_type: userextrapolation not implemented for spline interpolation",10);
      break;
    default:
      imthrow("translate_extrapolation_type: I am lost",10);
      break;
    }
    return(SPLINTERPOLATOR::Zeros);
  }

  // GENERAL MANIPULATION

  template <class T>
  void volume<T>::threshold(T lowerth, T upperth, threshtype tt)
  {
    if (activeROI) {
      for (int z=minz(); z<=maxz(); z++) {
	for (int y=miny(); y<=maxy(); y++) {
	  for (int x=minx(); x<=maxx(); x++) {
	    if ( ((tt==inclusive) && 
			(value(x,y,z)>= lowerth) && (value(x,y,z)<= upperth))  
	      || ((tt==exclusive) && 
			(value(x,y,z)>lowerth) && (value(x,y,z)<upperth)) ) 
	    {
	      //value(x,y,z) = 1;
	    } else {
	      value(x,y,z) = 0;
	    }
	  }
	}
      }
    } else {
      for (nonsafe_fast_iterator it=nsfbegin(), itend=nsfend();
	   it!=itend; ++it) {
	if ( ((tt==inclusive) && ((*it)>= lowerth) && ((*it)<= upperth)) 
	  || ((tt==exclusive) && ((*it)> lowerth) && ((*it)< upperth)) )
	{
	  //*it = 1;
	} else {
	  *it = 0;
	}
      }
    }
  }
  


  template <class T>
  void volume<T>::binarise(T lowerth, T upperth, threshtype tt)
  {
    if (activeROI) {
      for (int z=minz(); z<=maxz(); z++) {
	for (int y=miny(); y<=maxy(); y++) {
	  for (int x=minx(); x<=maxx(); x++) {
	    if ( ((tt==inclusive) && 
			(value(x,y,z)>= lowerth) && (value(x,y,z)<= upperth))  
	      || ((tt==exclusive) && 
			(value(x,y,z)>lowerth) && (value(x,y,z)<upperth)) ) 
	    {
	      value(x,y,z) = 1;
	    } else {
	      value(x,y,z) = 0;
	    }
	  }
	}
      }
    } else {
      for (nonsafe_fast_iterator it=nsfbegin(), itend=nsfend();
	   it!=itend; ++it) {
	if ( ((tt==inclusive) && ((*it)>= lowerth) && ((*it)<= upperth)) 
	  || ((tt==exclusive) && ((*it)> lowerth) && ((*it)< upperth)) )
	{
	  *it = 1;
	} else {
	  *it = 0;
	}
      }
    }
  }
  


  template <class S>
  inline S swapval(S xval, S yval, S zval, int dim)
  {
    switch (dim) {
    case 1:
      return xval;
    case -1:
      return -xval;
    case 2:
      return yval;
    case -2:
      return -yval;
    case 3:
      return zval;
    case -3:
      return -zval;
    }
    return (S) 0;  // should never get here
  }
 

  template <class T>
  inline int coordval(const volume<T>& vol, int x, int y, int z, int dim)
  {
    switch (dim) {
    case 1:
      return x;
    case -1:
      return vol.xsize() - x - 1;
    case 2:
      return y;
    case -2:
      return vol.ysize() - y - 1;
    case 3:
      return z;
    case -3:
      return vol.zsize() - z - 1;
    }
    return 0;  // should never get here
  }
  

  int dimarg(const string& val)
  {
    if (val=="x") {
      return 1;
    } else if (val=="x-" || val=="-x") {
      return -1;
    } else if (val=="y") {
      return 2;
    } else if (val=="y-" || val=="-y") {
      return -2;
    } else if (val=="z") {
      return 3;
    } else if (val=="z-" || val=="-z") {
      return -3;
    } else {
      return 0;
    }
  }
  

  template <class T>
  void setrow(Matrix& affmat, int rownum, int dimnum, const volume<T>& invol)
  {
    if (dimnum==1 || dimnum==-1) {
      affmat(rownum,1)=1*sign(dimnum); affmat(rownum,2)=0; affmat(rownum,3)=0;
    }
    if (dimnum==2 || dimnum==-2) {
      affmat(rownum,1)=0; affmat(rownum,2)=1*sign(dimnum); affmat(rownum,3)=0;
    }
    if (dimnum==3 || dimnum==-3) {
      affmat(rownum,1)=0; affmat(rownum,2)=0; affmat(rownum,3)=1*sign(dimnum);
    }
    if (dimnum>0) return;
    float fov=0.0;
    if (dimnum==-1) {
      fov = (invol.xsize() -1) * invol.xdim();
    }
    if (dimnum==-2) {
      fov = (invol.ysize() -1) * invol.ydim();
    }
    if (dimnum==-3) {
      fov = (invol.zsize() -1) * invol.zdim();
    }
    affmat(rownum,4)=fov;
  }
  

  template <class T>
  Matrix volume<T>::swapmat(int dim1, int dim2, int dim3) const
  {
    Matrix affmat(4,4);
    affmat = 0.0;
    affmat(4,4)=1.0;
    setrow(affmat,1,dim1,*this);
    setrow(affmat,2,dim2,*this);
    setrow(affmat,3,dim3,*this);
    return affmat;
  }


  template <class T>
  Matrix volume<T>::swapmat(const string& newx, const string& newy, const string& newz) const
  {
    return this->swapmat(dimarg(newx),dimarg(newy),dimarg(newz));
  }

  
  template <class T>
  void volume<T>::swapdimensions(const string& newx, const string& newy, const string& newz)
  {
    this->swapdimensions(dimarg(newx),dimarg(newy),dimarg(newz));
  }


  template <class T>
  void volume<T>::swapdimensions(int dim1, int dim2, int dim3)
  {
    basic_swapdimensions(dim1,dim2,dim3,true);
  }


  template <class T>
  void volume<T>::basic_swapdimensions(int dim1, int dim2, int dim3, bool keepLRorder)
  {
    // valid entries for dims are +/- 1, 2, 3 (corresponding to +/- x,y,z)
    if ( (dim1>3) || (dim1<-3) || (dim1==0) ||
	 (dim2<-3) || (dim2>3) || (dim2==0) || 
	 (dim3<-3) || (dim3>3) || (dim3==0) )
      { 
	imthrow("Invalid dimension numbers entered to swapdimensions",8);
      }

    if ( (std::abs(dim1)==std::abs(dim2)) || (std::abs(dim1)==std::abs(dim3)) 
	 || (std::abs(dim2)==std::abs(dim3)) ) 
      {
	imthrow("Dimension numbers were not a permutation in swapdimensions",8);
      }
    int sx, sy, sz;
    sx = std::abs(swapval(this->xsize(),this->ysize(),this->zsize(),dim1));
    sy = std::abs(swapval(this->xsize(),this->ysize(),this->zsize(),dim2));
    sz = std::abs(swapval(this->xsize(),this->ysize(),this->zsize(),dim3));
    volume<T> swapvol(sx,sy,sz);

    // now copy only the appropriate data
    for (int z=0; z<this->zsize(); z++) {
      for (int y=0; y<this->ysize(); y++) {
	for (int x=0; x<this->xsize(); x++) {
	  int nx = coordval(*this,x,y,z,dim1);
	  int ny = coordval(*this,x,y,z,dim2);
	  int nz = coordval(*this,x,y,z,dim3);
	  swapvol.value(nx,ny,nz) = this->value(x,y,z);
	}
      }
    }

    swapvol.copyproperties(*this);

    // now fix up the spatial properties

    // if a LR flip has happened then for all properties retain the unflipped values
    // therefore the data really has flipped, as otherwise it views identically
    if (keepLRorder) {
      // arbitrarily choose x to flip (if necessary)
      if (this->swapmat(dim1,dim2,dim3).Determinant() < 0) { dim1*=-1; }  
    }

    float dx = swapval(this->xdim(), this->ydim(), this->zdim(), dim1);
    float dy = swapval(this->xdim(), this->ydim(), this->zdim(), dim2);
    float dz = swapval(this->xdim(), this->ydim(), this->zdim(), dim3);
    swapvol.setdims(dx,dy,dz);

    // fix sform and qform matrices
    //   NB: sform and qform are voxel->mm but swapmat is mm->mm (flirt mm), 
    //   hence sampling mats

    Matrix nmat;
    nmat = this->sform_mat() * this->sampling_mat().i() * this->swapmat(dim1,dim2,dim3).i() * swapvol.sampling_mat();
    swapvol.set_sform(this->sform_code(), nmat);
    nmat = this->qform_mat() * this->sampling_mat().i() * this->swapmat(dim1,dim2,dim3).i() * swapvol.sampling_mat();
    swapvol.set_qform(this->qform_code(), nmat);
    
    int nx, ny, nz, mx, my, mz;
    mx = coordval(*this,this->minx(),this->miny(),this->minz(),dim1);
    my = coordval(*this,this->minx(),this->miny(),this->minz(),dim2);
    mz = coordval(*this,this->minx(),this->miny(),this->minz(),dim3);
    nx = coordval(*this,this->maxx(),this->maxy(),this->maxz(),dim1);
    ny = coordval(*this,this->maxx(),this->maxy(),this->maxz(),dim2);
    nz = coordval(*this,this->maxx(),this->maxy(),this->maxz(),dim3);
    swapvol.setROIlimits(mx,my,mz,nx,ny,nz);

    swapvol.set_whole_cache_validity(false);
    // now force the ROI and limits to be rebuilt
    swapvol.deactivateROI();
    if (this->usingROI())  swapvol.activateROI();

    *this = swapvol;
  }



  template <class T>
  int volume<T>::left_right_order() const
  {
    int order;
    // call the function in fslio
    order = FslGetLeftRightOrder2(this->sform_code(),
				   newmat_to_mat44(this->sform_mat()),
				   this->qform_code(),
				   newmat_to_mat44(this->qform_mat()));
    return order;
  }

  template <class T>
  short vox2mm_all(const volume<T>& vol, Matrix& vox2mm_mat, short& code) 
  {
    mat44 vox2mm44;
    code = FslGetVox2mmMatrix2(&vox2mm44,vol.sform_code(),
			      newmat_to_mat44(vol.sform_mat()), 
			      vol.qform_code(), 
			      newmat_to_mat44(vol.qform_mat()), 
			      vol.xdim(),vol.ydim(),vol.zdim());
    vox2mm_mat = mat44_to_newmat(vox2mm44);
    return code;
  }
 	  	 
  template <class T>
  Matrix volume<T>::newimagevox2mm_mat() const
  {
    Matrix vox2mm;
    short code;
    vox2mm_all(*this,vox2mm,code);
    return vox2mm;
  }

  template <class T>
  Matrix volume<T>::niftivox2newimagevox_mat() const
  {
    Matrix vox2vox=IdentityMatrix(4);
    if ((!RadiologicalFile) && (this->left_right_order()==FSL_RADIOLOGICAL)) {
      vox2vox = (this->sampling_mat()).i() * this->swapmat(-1,2,3) * 
	(this->sampling_mat());
    }
    return vox2vox;
  }
 	 


  template <class T>
  void volume<T>::swapLRorder()
  {
    basic_swapdimensions(-1,2,3,false);
  }


  template <class T>
  void volume<T>::setLRorder(int LRorder)
  {
    if (LRorder != this->left_right_order()) { this->swapLRorder(); }
  }


  template <class T>
  void volume<T>::makeradiological()
  {
    // use existing matrices to determine the order and if necessary swap
    //  all data and matrices
    if (this->left_right_order()==FSL_NEUROLOGICAL) { this->swapLRorder(); }
  }


  template <class T>
  void volume<T>::makeneurological()
  {
    // use existing matrices to determine the order and if necessary swap
    //  all data and matrices
    if (this->left_right_order()==FSL_RADIOLOGICAL) { this->swapLRorder(); }
  }

  


  // ARITHMETIC OPERATIONS


  template <class T>
  T volume<T>::operator=(T val)
  {
    if (activeROI) {
      for (int z=minz(); z<=maxz(); z++) {
	for (int y=miny(); y<=maxy(); y++) {
	  for (int x=minx(); x<=maxx(); x++) {
	    value(x,y,z) = val;
	  }
	}
      }
    } else {
      fill(nsfbegin(),nsfend(),val);  // use the STL
    }
    return val;
  }



  template <class T>
  const volume<T>& volume<T>::operator+=(T val)
  {
    if (activeROI) {
      for (int z=minz(); z<=maxz(); z++) {
	for (int y=miny(); y<=maxy(); y++) {
	  for (int x=minx(); x<=maxx(); x++) {
	    value(x,y,z) += val;
	  }
	}
      }
    } else {
      for (nonsafe_fast_iterator it=nsfbegin(), itend=nsfend();
	   it!=itend; ++it) {
	*it += val;
      }
    }
    return *this;
    
  }
    
  template <class T>
  const volume<T>& volume<T>::operator-=(T val)
  {
    if (activeROI) {
      for (int z=minz(); z<=maxz(); z++) {
	for (int y=miny(); y<=maxy(); y++) {
	  for (int x=minx(); x<=maxx(); x++) {
	    value(x,y,z) -= val;
	  }
	}
      }
    } else {
      for (nonsafe_fast_iterator it=nsfbegin(), itend=nsfend();
	   it!=itend; ++it) {
	*it -= val;
      }
    }
    return *this;
    
  }
    
  template <class T>
  const volume<T>& volume<T>::operator*=(T val)
  {
    if (activeROI) {
      for (int z=minz(); z<=maxz(); z++) {
	for (int y=miny(); y<=maxy(); y++) {
	  for (int x=minx(); x<=maxx(); x++) {
	    value(x,y,z) *= val;
	  }
	}
      }
    } else {
      for (nonsafe_fast_iterator it=nsfbegin(), itend=nsfend();
	   it!=itend; ++it) {
	*it *= val;
      }
    }
    return *this;
    
  }
    
  template <class T>
  const volume<T>& volume<T>::operator/=(T val)
  {
    if (activeROI) {
      for (int z=minz(); z<=maxz(); z++) {
	for (int y=miny(); y<=maxy(); y++) {
	  for (int x=minx(); x<=maxx(); x++) {
	    value(x,y,z) /= val;
	  }
	}
      }
    } else {
      for (nonsafe_fast_iterator it=nsfbegin(), itend=nsfend();
	   it!=itend; ++it) {
	*it /= val;
      }
    }
    return *this;
    
  }
    

  template <class T>
  const volume<T>& volume<T>::operator+=(const volume<T>& source)
  {
    if (!samesize(*this,source)) {
      imthrow("Attempted to add images/ROIs of different sizes",3);
    }
    if (activeROI || source.activeROI) {
      int xoff=source.minx()-minx(), yoff=source.miny()-miny(), 
	  zoff=source.minz()-minz();
      for (int z=minz(); z<=maxz(); z++) {
	for (int y=miny(); y<=maxy(); y++) {
	  for (int x=minx(); x<=maxx(); x++) {
	    value(x,y,z) += source(x+xoff,y+yoff,z+zoff);
	  }
	}
      }
    } else {
      fast_const_iterator dit=source.fbegin();
      for (nonsafe_fast_iterator it=nsfbegin(), itend=nsfend(); 
	   it!=itend; ++it, ++dit) {
	*it += *dit;
      }
    }
    return *this;
    
  }


  template <class T>
  const volume<T>& volume<T>::operator-=(const volume<T>& source)
  {
    if (!samesize(*this,source)) {
      imthrow("Attempted to subtract images/ROIs of different sizes",3);
    }
    if (activeROI || source.activeROI) {
      int xoff=source.minx()-minx(), yoff=source.miny()-miny(), 
	  zoff=source.minz()-minz();
      for (int z=minz(); z<=maxz(); z++) {
	for (int y=miny(); y<=maxy(); y++) {
	  for (int x=minx(); x<=maxx(); x++) {
	    value(x,y,z) -= source(x+xoff,y+yoff,z+zoff);
	  }
	}
      }
    } else {
      fast_const_iterator dit=source.fbegin();
      for (nonsafe_fast_iterator it=nsfbegin(), itend=nsfend(); 
	   it!=itend; ++it, ++dit) {
	*it -= *dit;
      }
    }
    return *this;
    
  }


  template <class T>
  const volume<T>& volume<T>::operator*=(const volume<T>& source)
  {
    if (!samesize(*this,source)) {
      imthrow("Attempted to multiply images/ROIs of different sizes",3);
    }
    if (activeROI || source.activeROI) {
      int xoff=source.minx()-minx(), yoff=source.miny()-miny(), 
	  zoff=source.minz()-minz();
      for (int z=minz(); z<=maxz(); z++) {
	for (int y=miny(); y<=maxy(); y++) {
	  for (int x=minx(); x<=maxx(); x++) {
	    value(x,y,z) *= source(x+xoff,y+yoff,z+zoff);
	  }
	}
      }
    } else {
      fast_const_iterator dit=source.fbegin();
      for (nonsafe_fast_iterator it=nsfbegin(), itend=nsfend(); 
	   it!=itend; ++it, ++dit) {
	*it *= *dit;
      }
    }
    return *this;
    
  }


  template <class T>
  const volume<T>& volume<T>::operator/=(const volume<T>& source)
  {
    if (!samesize(*this,source)) {
      imthrow("Attempted to divide images/ROIs of different sizes",3);
    }
    if (activeROI || source.activeROI) {
      int xoff=source.minx()-minx(), yoff=source.miny()-miny(), 
	  zoff=source.minz()-minz();
      for (int z=minz(); z<=maxz(); z++) {
	for (int y=miny(); y<=maxy(); y++) {
	  for (int x=minx(); x<=maxx(); x++) {
	    value(x,y,z) /= source(x+xoff,y+yoff,z+zoff);
	  }
	}
      }
    } else {
      fast_const_iterator dit=source.fbegin();
      for (nonsafe_fast_iterator it=nsfbegin(), itend=nsfend(); 
	   it!=itend; ++it, ++dit) {
	*it /= *dit;
      }
    }
    return *this;
    
  }


  template <class T>
  volume<T> volume<T>::operator+(T num) const
  {
    volume<T> tmp = *this;
    tmp+=num;
    return tmp;
  }

  template <class T>
  volume<T> volume<T>::operator-(T num) const
  {
    volume<T> tmp = *this;
    tmp-=num;
    return tmp;
  }

  template <class T>
  volume<T> volume<T>::operator*(T num) const
  {
    volume<T> tmp = *this;
    tmp*=num;
    return tmp;
  }

  template <class T>
  volume<T> volume<T>::operator/(T num) const
  {
    volume<T> tmp = *this;
    tmp/=num;
    return tmp;
  }

  template <class T>
  volume<T> volume<T>::operator+(const volume<T>& vol2) const
  {
    volume<T> tmp = *this;
    tmp+=vol2;
    return tmp;
  }

  template <class T>
  volume<T> volume<T>::operator-(const volume<T>& vol2) const
  {
    volume<T> tmp = *this;
    tmp-=vol2;
    return tmp;
  }

  template <class T>
  volume<T> volume<T>::operator*(const volume<T>& vol2) const
  {
    volume<T> tmp = *this;
    tmp*=vol2;
    return tmp;
  }

  template <class T>
  volume<T> volume<T>::operator/(const volume<T>& vol2) const
  {
    volume<T> tmp = *this;
    tmp/=vol2;
    return tmp;
  }




  ///////////////////////////////////////////////////////////////////////
  //////////////////////// VOLUME4D CLASS /////////////////////////////
  ///////////////////////////////////////////////////////////////////////

  // interfaces for lazy evaluated calculations
  template <class T>
  minmaxstuff<T> calc_minmax(const volume4D<T>& source);

  template <class T>
  minmaxstuff<T> calc_minmax(const volume4D<T>& source, const volume<T>& mask);

  template <class T>
  minmaxstuff<T> calc_minmax(const volume4D<T>& source, const volume4D<T>& mask);
 
  template <class T>
  std::vector<double> calc_sums(const volume4D<T>& vol);

  template <class T>
  std::vector<double> calc_sums(const volume4D<T>& vol, const volume<T>& mask);

  template <class T>
  std::vector<double> calc_sums(const volume4D<T>& vol, const volume4D<T>& mask);
 
  template <class T>
  std::vector<T> calc_percentiles(const volume4D<T>& vol);

  template <class T>
  std::vector<T> calc_percentiles(const volume4D<T>& vol, const volume<T>& mask, 
				  const std::vector<float>& percentilepvals);

  template <class T>
  std::vector<T> calc_percentiles(const volume4D<T>& vol, const volume4D<T>& mask, 
				  const std::vector<float>& percentilepvals);

  // CONSTRUCTORS (not including copy constructor - see under copying)

  template <class T>
  int volume4D<T>::initialize(int xsize, int ysize, int zsize, int tsize, 
				T *d)
  {
    this->destroy();
    long int vsize = xsize * ysize * zsize;
    volume<T> dummyvol;
    vols.insert(vols.begin(),tsize,dummyvol);
    T* ptr = d;
    for (int t=0; t<tsize; t++) {
      vols[t].reinitialize(xsize,ysize,zsize,ptr,false);
      if (ptr!=0)  ptr += vsize;
    }
    setdefaultproperties();
    return 0;
  }


  template <class T>
  void volume4D<T>::setdefaultproperties()
  {
    p_TR = 1.0;
    Limits.resize(8,0);
    setdefaultlimits();
    // Default ROI is whole volume
    ROIbox = Limits;
    activeROI = false;

    p_extrapmethod = zeropad;
    p_interpmethod = trilinear;
    p_padval = (T) 0;

    tsminmax.init(this,calc_minmax);
    sums.init(this,calc_sums);
    percentiles.init(this,calc_percentiles);
    robustlimits.init(this,calc_robustlimits);
    l_histogram.init(this,calc_histogram);

    // Initial percentile pvals to store when calculating percentiles
    percentilepvals.erase(percentilepvals.begin(),percentilepvals.end());
    percentilepvals.push_back(0.0);
    percentilepvals.push_back(0.001);
    percentilepvals.push_back(0.005);
    for (int probval=1; probval<=99; probval++) {
      percentilepvals.push_back(((float) probval)/100.0);
    }
    percentilepvals.push_back(0.995);
    percentilepvals.push_back(0.999);
    percentilepvals.push_back(1.0);    

    set_whole_cache_validity(false);
  }


  template <class T>
  volume4D<T>::volume4D()
    {
      this->initialize(0,0,0,0,0);
    }
  

  template <class T>
  volume4D<T>::volume4D(int xsize, int ysize, int zsize, int tsize, T *d) 
    {
      this->initialize(xsize,ysize,zsize,tsize,d);
    }


  template <class T>
  int volume4D<T>::reinitialize(int xsize, int ysize, int zsize, int tsize, 
				  T *d)
    {
      return this->initialize(xsize,ysize,zsize,tsize,d);
    }


  template <class T>
  void volume4D<T>::destroy()
    { 
      for (int t=0; t<tsize(); t++) { vols[t].destroy(); }
      if (tsize()>0)  vols.clear(); 
    }


  template <class T>
  volume4D<T>::~volume4D()
    { 
      this->destroy();
    }


  template <class T>
  void volume4D<T>::setdefaultlimits() const
    {
      Limits[0]=0; Limits[1]=0; Limits[2]=0; Limits[3]=0;
      Limits[4]=this->xsize()-1; 
      Limits[5]=this->ysize()-1; 
      Limits[6]=this->zsize()-1;
      Limits[7]=this->tsize()-1;
    }


  template <class T>
  void volume4D<T>::enforcelimits(std::vector<int>& lims) const
    {
//        lims[0]=Max(0,lims[0]); 
//        lims[1]=Max(0,lims[1]); 
//        lims[2]=Max(0,lims[2]); 
      lims[3]=Max(0,lims[3]); 
//        lims[4]=Min(this->xsize() - 1,lims[4]); 
//        lims[5]=Min(this->ysize() - 1,lims[5]); 
//        lims[6]=Min(this->zsize() - 1,lims[6]); 
      lims[7]=Min(this->tsize() - 1,lims[7]); 
    }

  // COPYING AND CONVERSION FUNCTIONS

  template <class T>
  int volume4D<T>::reinitialize(const volume4D<T>& source)
    {
      this->initialize(source.xsize(),source.ysize(),source.zsize(),
		       source.tsize(),0);
      this->copyvolumes(source);
      return this->copyproperties(source);
    }
  
  template <class T>
  volume4D<T>::volume4D(const volume4D<T>& source)
    {
      this->reinitialize(source);
    }

  template <class T>
  const volume4D<T>& volume4D<T>::operator=(const volume4D<T>& source)
  {
    this->reinitialize(source);
    return *this;
  }

  template <class T>
  const volume4D<T>& volume4D<T>::operator=(const volume<T>& source)
  {
    this->clear();
    this->addvolume(source);
    return *this;
  }

  template <class T>
  int volume4D<T>::copyvolumes(const volume4D<T>& source) {
    if (tsize() != source.tsize()) {
      imthrow("Attempted to copy with non-matching tsizes",2);
    }
    
    for (int t=0; t<source.tsize(); t++) {
      vols[t] = source.vols[t];
    }
    return 0;
  }


  template <class T>
  int volume4D<T>::copyproperties(const volume4D<T>& source) {
    // sets all properties
    copybasicproperties(source,*this);
    // copy lazy properties
    tsminmax.copy(source.tsminmax,this);
    sums.copy(source.sums,this);
    percentiles.copy(source.percentiles,this);
    percentilepvals = source.percentilepvals;
    robustlimits.copy(source.robustlimits,this);
    l_histogram.copy(source.l_histogram,this);
    HISTbins = source.HISTbins;
    HISTmin = source.HISTmin;
    HISTmax = source.HISTmax;

    // now copy all individual volume properties
    if (sameabssize(source,*this)) {
      for (int t=0; t<source.tsize(); t++) {
	vols[t].copyproperties(source[Min(source.tsize()-1,t)]);
      }
    } else {
      // Only copy ROIs
      int toffset = source.mint() - this->mint();
      for (int t=this->mint(); t<=this->maxt(); t++) {
	vols[t].copyproperties(source[Min(t + toffset,source.maxt())]);
      }
    }
    return 0;
  }


  template <class T>
  int volume4D<T>::copyproperties(const volume<T>& source) {
    // copy all individual volume properties
    for (int t=0; t<this->tsize(); t++) {
      vols[t].copyproperties(source);
    }
    return 0;
  }



  // ROI FUNCTIONS

  template <class T>
  volume4D<T> volume4D<T>::ROI() const
  {
    volume4D<T> roivol;
    roivol.reinitialize(this->maxx() - this->minx()+1,
			this->maxy() - this->miny()+1,
			this->maxz() - this->minz()+1,
			this->maxt() - this->mint()+1);

    // now copy only the appropriate data
    for (int t=this->mint(); t<=this->maxt(); t++) {
      roivol[t - this->mint()].copyROIonly(this->vols[t]);
    }
    roivol.copyproperties(*this);
    roivol.deactivateROI();
    // set sform and qform matrices appropriately (if set)
    Matrix roi2vol= IdentityMatrix(4);
    roi2vol(1,4) = this->minx();
    roi2vol(2,4) = this->miny();
    roi2vol(3,4) = this->minz();
    if (this->sform_code()!=NIFTI_XFORM_UNKNOWN) {
      roivol.set_sform(this->sform_code(),this->sform_mat() * roi2vol);
    }
    if (this->qform_code()!=NIFTI_XFORM_UNKNOWN) {
      roivol.set_qform(this->qform_code(),this->qform_mat() * roi2vol);
    }
    roivol.set_whole_cache_validity(false);
    return roivol;
  }


  template <class T>
  int volume4D<T>::copyROIonly(const volume4D<T>& source)
  {
    if (!samesize(*this,source)) {
      imthrow("Attempted to copy ROIs when different sizes",3);
    }

    int toff = mint() - source.mint();
    for (int t=source.mint(); t<=source.maxt(); t++) {
      vols[t + toff].copyROIonly(source[t]);
    }

    set_whole_cache_validity(false);  // is the default by reinitialize too
    return 0;
  }


  template <class T>
  void volume4D<T>::setROIlimits(int x0, int y0, int z0, int t0, 
				 int x1, int y1, int z1, int t1) const
  { 
    this->setROIlimits(x0,y0,z0,x1,y1,z1);
    this->setROIlimits(t0,t1);
  }
  

  template <class T>
  void volume4D<T>::setROIlimits(int x0, int y0, int z0, 
				 int x1, int y1, int z1) const
    { 
      // Enforce ordering 
      ROIbox[0]=Min(x0,x1); 
      ROIbox[1]=Min(y0,y1); 
      ROIbox[2]=Min(z0,z1); 
      ROIbox[4]=Max(x0,x1); 
      ROIbox[5]=Max(y0,y1); 
      ROIbox[6]=Max(z0,z1); 
      enforcelimits(ROIbox);
      for (int t=0; t<this->tsize(); t++) {
	vols[t].setROIlimits(x0,y0,z0,x1,y1,z1);
      }
      if (activeROI) this->activateROI();
    }

  template <class T>
  void volume4D<T>::setROIlimits(int t0, int t1) const
    { 
      // Enforce ordering 
      ROIbox[3]=Min(t0,t1); 
      ROIbox[7]=Max(t0,t1); 
      enforcelimits(ROIbox);
      if (activeROI) this->activateROI();
    }

  template <class T>
  void volume4D<T>::setROIlimits(const std::vector<int>& lims) const
      { 
	if (lims.size()!=8) return; 
	setROIlimits(lims[0],lims[1],lims[2],lims[3],
		     lims[4],lims[5],lims[6],lims[7]);
      }
 
  template <class T>
  void volume4D<T>::activateROI() const
  { 
    activeROI=true; 
    enforcelimits(ROIbox);
    Limits = ROIbox; 
    set_whole_cache_validity(false);
    for (int t=0; t<this->tsize(); t++) {
      vols[t].activateROI();
    }
  }
  
  template <class T>
  void volume4D<T>::deactivateROI() const
  { 
    activeROI=false; 
    setdefaultlimits();
    set_whole_cache_validity(false); 
    for (int t=0; t<this->tsize(); t++) {
      vols[t].deactivateROI();
    }
  }


  template <class T>
  void make_consistent_params(const volume4D<T>& vols, int t) 
    {
      vols[t].setextrapolationmethod(vols.getextrapolationmethod());
      vols[t].setinterpolationmethod(vols.getinterpolationmethod());
      if (vols.tsize()>0) vols[t].definekernelinterpolation(vols[0]);
      vols[t].setpadvalue(vols.getpadvalue());
      vols[t].setROIlimits(vols.minx(),vols.miny(),vols.minz(),
			   vols.maxx(),vols.maxy(),vols.maxz());
      if ( (vols[t].usingROI()) && (!vols.usingROI()) )
	{ vols[t].deactivateROI(); }
      if ( (!vols[t].usingROI()) && (vols.usingROI()) )
	{ vols[t].activateROI(); }
    }


  // VOLUME ACCESS FUNCTIONS (INSERT AND DELETE)
  
  template <class T>
  void volume4D<T>::insertvolume(const volume<T>& source, int t)
  {
    if ((t<0) || (t>tsize())) t=tsize();
    if (tsize()>0) {
      if (!sameabssize(source,vols[0])) {
	imthrow("Non-equal volume sizes in volume4D",3);
      }
    }
    vols.insert(vols.begin() + t,source);
    if (!activeROI) setdefaultlimits();  // fix up new time-length
    make_consistent_params(*this,t);
    set_whole_cache_validity(false);
  }

  template <class T>
  void volume4D<T>::addvolume(const volume<T>& source)
  { // add to the end of the series
    insertvolume(source,tsize()); 
  }

  template <class T>
  void volume4D<T>::addvolume(const volume4D<T>& source)
  { // add to the end of the series
    for (int t=source.mint(); t<=source.maxt(); t++) {
      addvolume(source[t]); 
    }
  }

  template <class T>
  void volume4D<T>::deletevolume(int t)
  {
    if ((t<0) || (t>tsize())) t=tsize();
    vols.erase(vols.begin() + t);
    if (!activeROI) setdefaultlimits();
    set_whole_cache_validity(false);
  }

  template <class T>
  void volume4D<T>::clear()
  {
    // inefficient but safe - should replace with a better method someday
    for (int t=tsize()-1; t>=0; t--) {
      this->deletevolume(t);
    }
  }



  template <class T>
  interpolation volume4D<T>::getinterpolationmethod() const
  { 
    return p_interpmethod;
  }

  template <class T>
  void volume4D<T>::setinterpolationmethod(interpolation interpmethod) const
  { 
    p_interpmethod = interpmethod;
    if (interpmethod == userinterpolation) {
      this->defineuserinterpolation(p_userinterp);
    }
    for (int t=0;  t<tsize(); t++) {
      vols[t].setinterpolationmethod(interpmethod);
      if ( (interpmethod == sinc) || (interpmethod == userkernel) ) {
	if (t>0) this->definekernelinterpolation(vols[0]);
      }
    }
  }

  template<class T>
  unsigned int volume4D<T>::getsplineorder() const
  {
    if (!tsize()) imthrow("getsplineorder: No volumes defined yet",10);
    return(vols[0].getsplineorder());
  }

  template<class T>
  void volume4D<T>::setsplineorder(unsigned int order) const
  {
    for (int i=0; i<tsize(); i++) vols[i].setsplineorder(order);
  }

  template<class T>
  void volume4D<T>::setextrapolationvalidity(bool xv, bool yv, bool zv) const
  {
    for (int i=0; i<tsize(); i++) vols[i].setextrapolationvalidity(xv,yv,zv);
  }

  template<class T>
  std::vector<bool> volume4D<T>::getextrapolationvalidity() const
  {
    if (!tsize()) imthrow("getextrapolationvalidity: No volumes defined yet",10);
    return(vols[0].getextrapolationvalidity());
  }  

  template <class T>
  extrapolation volume4D<T>::getextrapolationmethod() const
  { 
    return p_extrapmethod;
  }

  template <class T>
  void volume4D<T>::setextrapolationmethod(extrapolation extrapmethod) const
  { 
    p_extrapmethod = extrapmethod;
    for (int t=0;  t<tsize(); t++) vols[t].setextrapolationmethod(extrapmethod);
  }

  template <class T>
  void volume4D<T>::setpadvalue(T padval) const
  { 
    p_padval = padval;
    for (int t=0;  t<tsize(); t++) vols[t].setpadvalue(padval);
  }

  template <class T>
  T volume4D<T>::getpadvalue() const
  { 
    return p_padval;
  }

  template <class T>
  void volume4D<T>::defineuserinterpolation(float (*interp)(
                       const volume<T>& , float, float, float)) const
  { 
    p_userinterp = interp;
    for (int t=0;  t<tsize(); t++) vols[t].defineuserinterpolation(interp);
  }

  template <class T>
  void volume4D<T>::defineuserextrapolation(T (*extrap)(
                       const volume<T>& , int, int, int)) const
  { 
    p_userextrap = extrap;
    for (int t=0;  t<tsize(); t++) vols[t].defineuserextrapolation(extrap);
  }

  template <class T>
  void volume4D<T>::definekernelinterpolation(const ColumnVector& kx, 
					      const ColumnVector& ky,
					      const ColumnVector& kz, 
					      int wx, int wy, int wz) const
  { 
    // full widths
    for (int t=0;  t<tsize(); t++) 
      vols[t].definekernelinterpolation(kx,ky,kz,wx,wy,wz);

  }

  template <class T>
  void volume4D<T>::definekernelinterpolation(const volume4D<T>& vol) const
  { 
    if (vol.tsize()>0) {
      for (int t=0;  t<tsize(); t++) 
	vols[t].definekernelinterpolation(vol.vols[0]);
    }
  }

  template <class T>
  void volume4D<T>::definekernelinterpolation(const volume<T>& vol) const
  { 
    for (int t=0;  t<tsize(); t++) 
      vols[t].definekernelinterpolation(vols[0]);
  }

  template <class T>
  void volume4D<T>::definesincinterpolation(const string& sincwindowtype,
					    int w, int nstore) const
  { 
    // full width
    for (int t=0;  t<tsize(); t++) 
      vols[t].definesincinterpolation(sincwindowtype, w, nstore);
  }
  template <class T>
  void volume4D<T>::definesincinterpolation(const string& sincwindowtype,
					    int wx, int wy, int wz, 
					    int nstore) const
  { 
    // full widths
    for (int t=0;  t<tsize(); t++) 
      vols[t].definesincinterpolation(sincwindowtype, wx, wy, wz, nstore);
  }
  
  // MATRIX <-> VOLUME4D CONVERSIONS

  template <class T>
  ReturnMatrix volume4D<T>::matrix(const volume<T>& mask) const
  {
    Matrix matv;
    if (tsize()<=0) return matv;
    if (!samesize(mask,vols[0])) {
      imthrow("Mask of different size used in matrix()",3);
    }
    long cidx(1);
    matv.ReSize(this->maxt() - this->mint() + 1, no_mask_voxels(mask));
    int xoff = vols[0].minx() - mask.minx();
    int yoff = vols[0].miny() - mask.miny();
    int zoff = vols[0].minz() - mask.minz();
    int toff = 1 - this->mint();
    for (int z=mask.minz(); z<=mask.maxz(); z++) {
      for (int y=mask.miny(); y<=mask.maxy(); y++) {
	for (int x=mask.minx(); x<=mask.maxx(); x++) {
	  if (mask(x,y,z)>0) {
	    for (int t=this->mint(); t<=this->maxt(); t++) {
	      matv(t+toff,cidx) = vols[t](x+xoff,y+yoff,z+zoff);
	    }
	    cidx++;
	  }
	}
      }
    }
    matv.Release();
    return matv;
  }

  template <class T>
  ReturnMatrix volume4D<T>::matrix(const volume<T>& mask, vector<long>& voxelLabels) const
  {
    voxelLabels.clear();
    Matrix matv;
    if (tsize()<=0) return matv;
    if (!samesize(mask,vols[0])) {
      imthrow("Mask of different size used in matrix()",3);
    }
    long cidx (1);
    matv.ReSize(this->maxt() - this->mint() + 1, no_mask_voxels(mask) );
    int xoff = vols[0].minx() - mask.minx();
    int yoff = vols[0].miny() - mask.miny();
    int zoff = vols[0].minz() - mask.minz();
    int toff = 1 - this->mint();
    for (int z=mask.minz(); z<=mask.maxz(); z++) {
      for (int y=mask.miny(); y<=mask.maxy(); y++) {
	for (int x=mask.minx(); x<=mask.maxx(); x++) {
	  if (mask(x,y,z)>0) {
	    voxelLabels.push_back(x+y*mask.xsize()+z*mask.xsize()*mask.ysize());
	    for (int t=this->mint(); t<=this->maxt(); t++) {
	      matv(t+toff,cidx) = vols[t](x+xoff,y+yoff,z+zoff);
	    }
	    cidx++;
	  }
	}
      }
    }
    matv.Release();
    return matv;
  }


  template <class T>
  ReturnMatrix volume4D<T>::matrix() const
  {
    volume<T> mask(vols[0]);
    mask = 1;
    return matrix(mask);
  }


  template <class T>
  void volume4D<T>::setmatrix(const Matrix& newmatrix, const volume<T>& mask,
			      const T pad)
  {
    int tsz = this->maxt() - this->mint() + 1;
    if ( (tsz==0) || 
	 (tsz!=newmatrix.Nrows()) || (!samesize(mask,vols[0])) ) {
      this->reinitialize(mask.xsize(),mask.ysize(),mask.zsize(),
			 newmatrix.Nrows());
    }
    this->copyproperties(mask);		       
    this->operator=(pad);
    if (newmatrix.Ncols()!=no_mask_voxels(mask)) {
      imthrow("Incompatible number of mask positions and matrix columns",4);
    }
    long cidx = 1;
    int xoff = mask.minx() - vols[0].minx();
    int yoff = mask.miny() - vols[0].miny();
    int zoff = mask.minz() - vols[0].minz();
    for (int z=vols[0].minz(); z<=vols[0].maxz(); z++) {
      for (int y=vols[0].miny(); y<=vols[0].maxy(); y++) {
	for (int x=vols[0].minx(); x<=vols[0].maxx(); x++) {
	  if (mask(x+xoff,y+yoff,z+zoff)>0) {
	    for (int t=this->mint(); t<=this->maxt(); t++) {
	      vols[t](x,y,z) = (T) newmatrix(t+1,cidx);
	    }
	    cidx++;
	  }
	}
      }
    }

    set_whole_cache_validity(false);
  }
      

  template <class T>
  void volume4D<T>::setmatrix(const Matrix& newmatrix)
  {
    volume<T> dummymask(vols[0]);
    dummymask = 1;
    this->setmatrix(newmatrix,dummymask,0);
  }
  
  template <class T>
  volume<int> volume4D<T>::vol2matrixkey(volume<T>& mask)
  {
    int count=1;
    volume<int> tmp(this->xsize(),this->ysize(),this->zsize());
    for(int z=0;z< this->zsize();z++){
      for(int y=0;y<this->ysize();y++){
	for(int x=0;x<this->xsize();x++){
	  if(mask(x,y,z)>0){
	    tmp(x,y,z)=count;
	    count++;
	  }
	  else{
	    tmp(x,y,z)=0;
	  }
	  
	}
      }
    }
   
    return tmp;
  }

  template <class T>
  ReturnMatrix volume4D<T>::matrix2volkey(volume<T>& mask){
    int count=0;
    for(int z=0;z< this->zsize();z++)
      for(int y=0;y<this->ysize();y++)
	for(int x=0;x<this->xsize();x++)
	  if(mask(x,y,z)>0)
	    count++;

    Matrix key(count,3);
    count=1;
    for(int z=0;z< this->zsize();z++)
      for(int y=0;y<this->ysize();y++)
	for(int x=0;x<this->xsize();x++)
	  if(mask(x,y,z)>0){
	    key(count,1)=x;
	    key(count,2)=y;
	    key(count,3)=z;
	    count++;
	  }
    key.Release();
    return key;
     
  }




  template <class T>
  ReturnMatrix volume4D<T>::voxelts(int x, int y, int z) const
  {
    ColumnVector res;
    if (this->maxt()<this->mint()) return res;
    res.ReSize(this->maxt() - this->mint() + 1);
    int toff = 1 - this->mint();
    for (int t=this->mint(); t<=this->maxt(); t++) {
      res(t + toff) = (NEWMAT::Real) vols[t](x,y,z);
    }

    res.Release();
    return res;
  }

  template <class T>
  void volume4D<T>::setvoxelts(const ColumnVector& ts, int x, int y, int z)
  {
    if (ts.Nrows() != (this->maxt() - this->mint() + 1)) {
      imthrow("setvoxelts - incorrectly sized vector",3);
    }
    for (int t=this->mint(); t<=this->maxt(); t++) {
      vols[t](x,y,z) = (T) ts(t+1);
    }
  }



  // PROPERTIES

  template <class T>
  void volume4D<T>::setxdim(float x) 
    { for (int t=0; t<tsize(); t++) vols[t].setxdim(x); }

  template <class T>
  void volume4D<T>::setydim(float y)
    { for (int t=0; t<tsize(); t++) vols[t].setydim(y); }

  template <class T>
  void volume4D<T>::setzdim(float z)
    { for (int t=0; t<tsize(); t++) vols[t].setzdim(z); }

  template <class T>
  void volume4D<T>::setDisplayMaximumMinimum(const float maximum, const float minimum) const
  { for (int t=0; t<tsize(); t++) vols[t].setDisplayMaximumMinimum(maximum,minimum); }

  template <class T>
  void volume4D<T>::setAuxFile(const string fileName)
  { for (int t=0; t<tsize(); t++) vols[t].setAuxFile(fileName); }

  // SECONDARY PROPERTIES (using a 3D or 4D mask)

  template <class T>
  minmaxstuff<T> calc_minmax(const volume4D<T>& source)
  {
    minmaxstuff<T> newminmax;
    newminmax.min=source(source.minx(),source.miny(),source.minz(),source.mint());
    newminmax.max=newminmax.min;
    newminmax.minx=source.minx();
    newminmax.miny=source.miny(); 
    newminmax.minz=source.minz(); 
    newminmax.maxx=source.minx(); 
    newminmax.maxy=source.miny(); 
    newminmax.maxz=source.minz();
    newminmax.mint=source.mint(); 
    newminmax.maxt=source.maxt();


    if (source.maxt()>=source.mint()) {
      newminmax = calc_minmax(source[source.mint()]);
      newminmax.mint = source.mint();
      newminmax.maxt = source.mint();
    }
    for (int t=source.mint(); t<=source.maxt(); t++) {
      if (source[t].min() < newminmax.min) { 
	newminmax.min = source[t].min();
	newminmax.minx = source[t].mincoordx();
	newminmax.miny = source[t].mincoordy();
	newminmax.minz = source[t].mincoordz();
	newminmax.mint = t;
      }
      if (source[t].max() > newminmax.max) {
	newminmax.max = source[t].max();
	newminmax.maxx = source[t].maxcoordx();
	newminmax.maxy = source[t].maxcoordy();
	newminmax.maxz = source[t].maxcoordz();
	newminmax.maxt = t;
      }
    }
    return newminmax;
  }


  template <class T>
  minmaxstuff<T> calc_minmax(const volume4D<T>& source, const volume<T>& mask)
  {
    if (!samesize(source[0],mask)) {
      imthrow("Mask of different size used in calc_minmax",3);
    }
    
    minmaxstuff<T> newminmax;
    newminmax.min=source(source.minx(),source.miny(),source.minz(),source.mint());
    newminmax.max=newminmax.min;
    newminmax.minx=source.minx();
    newminmax.miny=source.miny(); 
    newminmax.minz=source.minz(); 
    newminmax.maxx=source.minx(); 
    newminmax.maxy=source.miny(); 
    newminmax.maxz=source.minz();
    newminmax.mint=source.mint(); 
    newminmax.maxt=source.maxt();
    
    
    if (source.maxt()>=source.mint()) {
      newminmax = calc_minmax(source[source.mint()],mask);
      newminmax.mint = source.mint();
      newminmax.maxt = source.mint();
    }
    for (int t=source.mint(); t<=source.maxt(); t++) {
      if (source[t].min(mask) < newminmax.min) { 
	newminmax.min = source[t].min(mask);
	newminmax.minx = source[t].mincoordx(mask);
	newminmax.miny = source[t].mincoordy(mask);
	newminmax.minz = source[t].mincoordz(mask);
	newminmax.mint = t;
      }
      if (source[t].max(mask) > newminmax.max) {
	newminmax.max = source[t].max(mask);
	newminmax.maxx = source[t].maxcoordx(mask);
	newminmax.maxy = source[t].maxcoordy(mask);
	newminmax.maxz = source[t].maxcoordz(mask);
	newminmax.maxt = t;
      }
    }
    return newminmax;
  }
  
  

  template <class T>
  minmaxstuff<T> calc_minmax(const volume4D<T>& source, const volume4D<T>& mask)
  {
    if (!samesize(source[0],mask[0])) {
      imthrow("Mask of different size used in calc_minmax",3);
    }
    
    minmaxstuff<T> newminmax;
    newminmax.min=source(source.minx(),source.miny(),source.minz(),source.mint());
    newminmax.max=newminmax.min;
    newminmax.minx=source.minx();
    newminmax.miny=source.miny(); 
    newminmax.minz=source.minz(); 
    newminmax.maxx=source.minx(); 
    newminmax.maxy=source.miny(); 
    newminmax.maxz=source.minz();
    newminmax.mint=source.mint(); 
    newminmax.maxt=source.maxt();


    if (source.maxt()>=source.mint()) {
      newminmax = calc_minmax(source[source.mint()],mask[Min(source.mint(),mask.maxt())]);
      newminmax.mint = source.mint();
      newminmax.maxt = source.mint();
    }
    for (int t=source.mint(); t<=source.maxt(); t++) {
      if (source[t].min(mask[Min(t,mask.maxt())]) < newminmax.min) { 
	newminmax.min = source[t].min(mask[Min(t,mask.maxt())]);
	newminmax.minx = source[t].mincoordx(mask[Min(t,mask.maxt())]);
	newminmax.miny = source[t].mincoordy(mask[Min(t,mask.maxt())]);
	newminmax.minz = source[t].mincoordz(mask[Min(t,mask.maxt())]);
	newminmax.mint = t;
      }
      if (source[t].max(mask[Min(t,mask.maxt())]) > newminmax.max) {
	newminmax.max = source[t].max(mask[Min(t,mask.maxt())]);
	newminmax.maxx = source[t].maxcoordx(mask[Min(t,mask.maxt())]);
	newminmax.maxy = source[t].maxcoordy(mask[Min(t,mask.maxt())]);
	newminmax.maxz = source[t].maxcoordz(mask[Min(t,mask.maxt())]);
	newminmax.maxt = t;
      }
    }
    return newminmax;
  }
  
  
  template <class T>
  std::vector<double> calc_sums(const volume4D<T>& vol)
  {
    std::vector<double> newsums(2), addterm(2);
    newsums[0]=0;  newsums[1]=0;
    for (int t=vol.mint(); t<=vol.maxt(); t++) {
      addterm = calc_sums(vol[t]);
      newsums[0] += addterm[0];
      newsums[1] += addterm[1];
    }
    return newsums;
  }


  template <class T>
  std::vector<double> calc_sums(const volume4D<T>& vol, const volume<T>& mask)
  {
    if (!samesize(vol[0],mask)) {
      imthrow("calc_sums:: mask and volume must be the same size",4);
    }
    std::vector<double> newsums(2), addterm(2);
    newsums[0]=0;  newsums[1]=0;
    for (int t=vol.mint(); t<=vol.maxt(); t++) {
      addterm = calc_sums(vol[t],mask);
      newsums[0] += addterm[0];
      newsums[1] += addterm[1];
    }
    return newsums;
  }

  
  template <class T>
  std::vector<double> calc_sums(const volume4D<T>& vol, const volume4D<T>& mask)
  {
    if (!samesize(vol[0],mask[0])) {
      imthrow("calc_sums:: mask and volume must be the same size",4);
    }
    std::vector<double> newsums(2), addterm(2);
    newsums[0]=0;  newsums[1]=0;
    for (int t=vol.mint(); t<=vol.maxt(); t++) {
      addterm = calc_sums(vol[t],mask[Min(t,mask.maxt())]);
      newsums[0] += addterm[0];
      newsums[1] += addterm[1];
    }
    return newsums;
  }


  
  template <class T>
  T volume4D<T>::min(const volume<T>& mask) const
  {
    minmaxstuff<T> retval;
    retval = calc_minmax(*this,mask);
    return retval.min;	
  }
  
  template <class T>
  T volume4D<T>::min(const volume4D<T>& mask) const
  {
    minmaxstuff<T> retval;
    retval = calc_minmax(*this,mask);
    return retval.min;	
  }

  template <class T>
  T volume4D<T>::max(const volume<T>& mask) const
  {
    minmaxstuff<T> retval;
    retval = calc_minmax(*this,mask);
    return retval.max;	
  }
  
  template <class T>
  T volume4D<T>::max(const volume4D<T>& mask) const
  {
    minmaxstuff<T> retval;
    retval = calc_minmax(*this,mask);
    return retval.max;	
  }

  template <class T>
  int volume4D<T>::mincoordx(const volume<T>& mask) const
  {
    minmaxstuff<T> retval;
    retval = calc_minmax(*this,mask);
    return retval.minx;	
  }
  
  template <class T>
  int volume4D<T>::mincoordx(const volume4D<T>& mask) const
  {
    minmaxstuff<T> retval;
    retval = calc_minmax(*this,mask);
    return retval.minx;	
  }
  
  template <class T>
  int volume4D<T>::mincoordy(const volume<T>& mask) const
  {
    minmaxstuff<T> retval;
    retval = calc_minmax(*this,mask);
    return retval.miny;	
  }
  
  template <class T>
  int volume4D<T>::mincoordy(const volume4D<T>& mask) const
  {
    minmaxstuff<T> retval;
    retval = calc_minmax(*this,mask);
    return retval.miny;	
  }
  
  template <class T>
  int volume4D<T>::mincoordz(const volume<T>& mask) const
  {
    minmaxstuff<T> retval;
    retval = calc_minmax(*this,mask);
    return retval.minz;	
  }
  
  template <class T>
  int volume4D<T>::mincoordz(const volume4D<T>& mask) const
  {
    minmaxstuff<T> retval;
    retval = calc_minmax(*this,mask);
    return retval.minz;	
  }
  
  template <class T>
  int volume4D<T>::maxcoordx(const volume<T>& mask) const
  {
    minmaxstuff<T> retval;
    retval = calc_minmax(*this,mask);
    return retval.maxx;	
  }
  
  template <class T>
  int volume4D<T>::maxcoordx(const volume4D<T>& mask) const
  {
    minmaxstuff<T> retval;
    retval = calc_minmax(*this,mask);
    return retval.maxx;	
  }
  
  template <class T>
  int volume4D<T>::maxcoordy(const volume<T>& mask) const
  {
    minmaxstuff<T> retval;
    retval = calc_minmax(*this,mask);
    return retval.maxy;	
  }
  
  template <class T>
  int volume4D<T>::maxcoordy(const volume4D<T>& mask) const
  {
    minmaxstuff<T> retval;
    retval = calc_minmax(*this,mask);
    return retval.maxy;	
  }
  
  template <class T>
  int volume4D<T>::maxcoordz(const volume<T>& mask) const
  {
    minmaxstuff<T> retval;
    retval = calc_minmax(*this,mask);
    return retval.maxz;	
  }
  
  template <class T>
  int volume4D<T>::maxcoordz(const volume4D<T>& mask) const
  {
    minmaxstuff<T> retval;
    retval = calc_minmax(*this,mask);
    return retval.maxz;	
  }

  template <class T>
  double volume4D<T>::sum(const volume<T>& mask) const
  {
    std::vector<double> retval;
    retval = calc_sums(*this,mask);
    return retval[0];	
  }

  template <class T>
  double volume4D<T>::sum(const volume4D<T>& mask) const
  {
    std::vector<double> retval;
    retval = calc_sums(*this,mask);
    return retval[0];	
  }

  template <class T>
  double volume4D<T>::sumsquares(const volume<T>& mask) const
  {
    std::vector<double> retval;
    retval = calc_sums(*this,mask);
    return retval[1];	
  }

  template <class T>
  double volume4D<T>::sumsquares(const volume4D<T>& mask) const
  {
    std::vector<double> retval;
    retval = calc_sums(*this,mask);
    return retval[1];	
  }

  template <class T>
  double volume4D<T>::mean(const volume<T>& mask) const
  { 
    return sum(mask)/(Max((double) no_mask_voxels(mask),1.0));
  }

  template <class T>
  double volume4D<T>::mean(const volume4D<T>& mask) const
  { 
    return sum(mask)/(Max((double) no_mask_voxels(mask),1.0));
  }


  template <class T>
  double volume4D<T>::variance(const volume<T>& mask) const
  { 
    if (no_mask_voxels(mask)>0) {
      double n=(double) no_mask_voxels(mask);
      return (n/Max(1.0,n-1))*(sumsquares(mask)/n - mean(mask)*mean(mask));
    } else {
      cerr << "ERROR:: Empty mask image" << endl;
      return 0;
    }
  }

  template <class T>
  double volume4D<T>::variance(const volume4D<T>& mask) const
  { 
    if (no_mask_voxels(mask)>0) {
      double n=(double) no_mask_voxels(mask);
      return (n/Max(1.0,n-1))*(sumsquares(mask)/n - mean(mask)*mean(mask));
    } else {
      cerr << "ERROR:: Empty mask image" << endl;
      return 0;
    }
  }



  template <class T>
  T volume4D<T>::percentile(float pvalue) const
  {
    if ((pvalue>1.0) || (pvalue<0.0)) 
      { imthrow("Percentiles must be in the range [0.0,1.0]",4); }
    int idx = get_pval_index(percentilepvals,pvalue);
    if (idx==pval_index_end()) {
      percentilepvals.push_back(pvalue);
      idx = percentilepvals.size() - 1;
      percentiles.force_recalculation();
    }
    assert((idx>=0) && (idx < (int) percentilepvals.size()));
    return percentiles()[idx];
  }


  template <class T>
  T volume4D<T>::percentile(float pvalue, const volume<T>& mask) const
  {
    if ((pvalue>1.0) || (pvalue<0.0)) 
      { imthrow("Percentiles must be in the range [0.0,1.0]",4); }
    std::vector<float> pvaluevec;
    std::vector<T> retval;
    pvaluevec.push_back(pvalue);
    retval = calc_percentiles(*this,mask,pvaluevec);
    return retval[0];
  }


  template <class T>
  T volume4D<T>::percentile(float pvalue, const volume4D<T>& mask) const
  {
    if ((pvalue>1.0) || (pvalue<0.0)) 
      { imthrow("Percentiles must be in the range [0.0,1.0]",4); }
    std::vector<float> pvaluevec;
    std::vector<T> retval;
    pvaluevec.push_back(pvalue);
    retval = calc_percentiles(*this,mask,pvaluevec);
    return retval[0];
  }



  template <class T>
  std::vector<T> calc_percentiles(const volume4D<T>& vol)
  {
    unsigned int numbins = (unsigned int) vol.nvoxels() * vol.ntimepoints();
    unsigned int hindx = 0;
    std::vector<T> hist(numbins);
    for (int t=vol.mint(); t<=vol.maxt(); t++) {
      for (int z=vol.minz(); z<=vol.maxz(); z++) {
	for (int y=vol.miny(); y<=vol.maxy(); y++) {
	  for (int x=vol.minx(); x<=vol.maxx(); x++) {
	    hist[hindx++] = vol(x,y,z,t);
	  }
	}
      }
    }
    return percentile_vec(hist,vol.percentilepvalues());
  }      

  
  template <class T>
  std::vector<T> calc_percentiles(const volume4D<T>& vol, const volume<T>& mask, 
				    const std::vector<float>& percentilepvals)
  {
    if (!samesize(vol[0],mask)) {
      imthrow("mask and vol have different sizes in calc_percentiles",3);
    }
    std::vector<T> hist;
    for (int z=vol.minz(); z<=vol.maxz(); z++) {
      for (int y=vol.miny(); y<=vol.maxy(); y++) {
	for (int x=vol.minx(); x<=vol.maxx(); x++) {
	  if (mask(x,y,z)>0.5) {
	    for (int t=vol.mint(); t<=vol.maxt(); t++) {
	      hist.push_back(vol(x,y,z,t));
	    }
	  }
	}
      }
    }
    return percentile_vec(hist,percentilepvals);
  }      
  

  template <class T>
  std::vector<T> calc_percentiles(const volume4D<T>& vol, const volume4D<T>& mask, 
				    const std::vector<float>& percentilepvals)
  {
    if (!samesize(vol[0],mask[0])) {
      imthrow("mask and vol have different sizes in calc_percentiles",3);
    }
    std::vector<T> hist;
    for (int t=vol.mint(); t<=vol.maxt(); t++) {
      for (int z=vol.minz(); z<=vol.maxz(); z++) {
	for (int y=vol.miny(); y<=vol.maxy(); y++) {
	  for (int x=vol.minx(); x<=vol.maxx(); x++) {
	    if (mask(x,y,z,Min(t,mask.maxt()))>0.5) hist.push_back(vol(x,y,z,t));
	  }
	}
      }
    }
    return percentile_vec(hist,percentilepvals);
  }      



  template <class T>
  ColumnVector calc_histogram(const volume4D<T>& vol)
  {
    ColumnVector hist;
    calc_histogram(vol,vol.histbins(),vol.histmin(),vol.histmax(),hist);
    return hist;
  }


  template <class T>
  ColumnVector calc_histogram(const volume4D<T>& vol, const volume4D<T>& mask)
  {
    ColumnVector hist;
    calc_histogram(vol,vol.histbins(),vol.histmin(),vol.histmax(),hist,mask);
    return hist;
  }


  template <class T>
  ColumnVector calc_histogram(const volume4D<T>& vol, const volume<T>& mask)
  {
    ColumnVector hist;
    calc_histogram(vol,vol.histbins(),vol.histmin(),vol.histmax(),hist,mask);
    return hist;
  }

 
  template <class T>
  ColumnVector volume4D<T>::histogram(int nbins, T minval, T maxval) const
  {
    bool sameparams = true;
    if (HISTbins != nbins) {
      HISTbins = nbins;
      sameparams = false;
    }
    if (HISTmin != minval) {
      HISTmin = minval;
      sameparams = false;
    }
    if (HISTmax != maxval) {
      HISTmax = maxval;
      sameparams = false;
    }
    if (!sameparams) {
      l_histogram.force_recalculation();
    }
    return l_histogram();
  }


  template <class T>
  ColumnVector volume4D<T>::histogram(int nbins) const
  {
    return histogram(nbins,robustmin(),robustmax());
  }

  template <class T>
  ColumnVector volume4D<T>::histogram(int nbins, T minval, T maxval, 
				      const volume4D<T>& mask) const
  {
    ColumnVector hist;
    calc_histogram(*this,nbins,minval,maxval,hist,mask);
    return hist;
  }

  template <class T>
  ColumnVector volume4D<T>::histogram(int nbins, const volume4D<T>& mask) const
  {
    return histogram(nbins,robustmin(),robustmax(),mask);
  }
    
  template <class T>
  ColumnVector volume4D<T>::histogram(int nbins, T minval, T maxval, 
				      const volume<T>& mask) const
  {
    ColumnVector hist;
    calc_histogram(*this,nbins,minval,maxval,hist,mask);
    return hist;
  }

  template <class T>
  ColumnVector volume4D<T>::histogram(int nbins, const volume<T>& mask) const
  {
    return histogram(nbins,robustmin(),robustmax(),mask);
  }
    
  

 // GENERAL MANIPULATION

  template <class T>
  void volume4D<T>::binarise(T lowerth, T upperth, threshtype tt)
  {
    set_whole_cache_validity(false);
    for (int t=this->mint(); t<=this->maxt(); t++)  {
      vols[t].binarise(lowerth,upperth,tt); 
    }
  }
  
  template <class T>
  void volume4D<T>::threshold(T lowerth, T upperth, threshtype tt)
  {
    set_whole_cache_validity(false);
    for (int t=this->mint(); t<=this->maxt(); t++)  { 
      vols[t].threshold(lowerth,upperth,tt); 
    }
  }
  
  template <class T>
  void volume4D<T>::swapdimensions(int dim1, int dim2, int dim3)
  {
    for (int t=0; t<this->tsize(); t++) {
      vols[t].swapdimensions(dim1,dim2,dim3);
    }
  }

  template <class T>
  void volume4D<T>::swapdimensions(const string& newx, const string& newy, const string& newz)
  {
    this->swapdimensions(dimarg(newx),dimarg(newy),dimarg(newz));
  }


  template <class T>
  Matrix volume4D<T>::swapmat(int dim1, int dim2, int dim3) const
  {
    if (this->tsize()>0) {
      return vols[0].swapmat(dim1,dim2,dim3);
    }
    return IdentityMatrix(4);
  }


  template <class T>
  Matrix volume4D<T>::swapmat(const string& newx, const string& newy, const string& newz) const
  {
    return this->swapmat(dimarg(newx),dimarg(newy),dimarg(newz));
  }


  template <class T>
  Matrix volume4D<T>::newimagevox2mm_mat() const
  {
    if (this->tsize()>0) return vols[0].newimagevox2mm_mat();
    return IdentityMatrix(4);
  }
 	 
  template <class T>
  Matrix volume4D<T>::niftivox2newimagevox_mat() const
  {
    if (this->tsize()>0) return vols[0].niftivox2newimagevox_mat();
    return IdentityMatrix(4);
  }

 	 
  template <class T>
  int volume4D<T>::left_right_order() const
  {
    if (this->tsize()>0) return vols[0].left_right_order();
    return FSL_RADIOLOGICAL;
  }

  template <class T>
  void volume4D<T>::swapLRorder()
  {
    for (int t=0; t<this->tsize(); t++) {
      vols[t].swapLRorder();
    }
  }


  template <class T>
  void volume4D<T>::setLRorder(int LRorder)
  {
    if (LRorder != this->left_right_order()) { this->swapLRorder(); }
  }


  template <class T>
  void volume4D<T>::makeradiological()
  {
    if (this->left_right_order()==FSL_NEUROLOGICAL) { this->swapLRorder(); }
  }


  template <class T>
  void volume4D<T>::makeneurological()
  {
    if (this->left_right_order()==FSL_RADIOLOGICAL) { this->swapLRorder(); }
  }

  // ARITHMETIC OPERATIONS


  template <class T>
  T volume4D<T>::operator=(T val)
  {
    set_whole_cache_validity(false);
    for (int t=this->mint(); t<=this->maxt(); t++)  { vols[t] = val; }
    return val;
  }


  template <class T>
  const volume4D<T>& volume4D<T>::operator+=(T val)
  {
    set_whole_cache_validity(false);
    for (int t=this->mint(); t<=this->maxt(); t++) {
      vols[t] += val;
    }
    return *this;
  }
    
  template <class T>
  const volume4D<T>& volume4D<T>::operator-=(T val)
  {
    set_whole_cache_validity(false);
    for (int t=this->mint(); t<=this->maxt(); t++) {
      vols[t] -= val;
    }
    return *this;
  }
    
  template <class T>
  const volume4D<T>& volume4D<T>::operator*=(T val)
  {
    set_whole_cache_validity(false);
    for (int t=this->mint(); t<=this->maxt(); t++) {
      vols[t] *= val;
    }
    return *this;
  }
    
  template <class T>
  const volume4D<T>& volume4D<T>::operator/=(T val)
  {
    set_whole_cache_validity(false);
    for (int t=this->mint(); t<=this->maxt(); t++) {
      vols[t] /= val;
    }
    return *this;
  }
    

  template <class T>
  const volume4D<T>& volume4D<T>::operator+=(const volume<T>& source)
  {
    set_whole_cache_validity(false);
    for (int t=this->mint(); t<=this->maxt(); t++) {
      vols[t] += source;
    }
    return *this;
  }

  template <class T>
  const volume4D<T>& volume4D<T>::operator-=(const volume<T>& source) 
  {
    set_whole_cache_validity(false);
    for (int t=this->mint(); t<=this->maxt(); t++) {
      vols[t] -= source;
    }
    return *this;
  }

  template <class T>
  const volume4D<T>& volume4D<T>::operator*=(const volume<T>& source) 
  {
    set_whole_cache_validity(false);
    for (int t=this->mint(); t<=this->maxt(); t++) {
      vols[t] *= source;
    }
    return *this;
  }

  template <class T>
  const volume4D<T>& volume4D<T>::operator/=(const volume<T>& source) 
  {
    set_whole_cache_validity(false);
    for (int t=this->mint(); t<=this->maxt(); t++) {
      vols[t] /= source;
    }
    return *this;
  }

 
  template <class T>
  const volume4D<T>& volume4D<T>::operator+=(const volume4D<T>& source)
  {
    if (!samesize(*this,source)) {
      imthrow("Attempted to add images/ROIs of different sizes",3);
    }
    set_whole_cache_validity(false);
    int toff = source.mint() - this->mint();
    for (int t=this->mint(); t<=this->maxt(); t++) {
      vols[t] += source[t + toff];
    }
    return *this;
  }


  template <class T>
  const volume4D<T>& volume4D<T>::operator-=(const volume4D<T>& source)
  {
    if (!samesize(*this,source)) {
      imthrow("Attempted to subtract images/ROIs of different sizes",3);
    }
    set_whole_cache_validity(false);
    int toff = source.mint() - this->mint();
    for (int t=this->mint(); t<=this->maxt(); t++) {
      vols[t] -= source[t + toff];
    }
    return *this;
  }


  template <class T>
  const volume4D<T>& volume4D<T>::operator*=(const volume4D<T>& source)
  {
    if (!samesize(*this,source)) {
      imthrow("Attempted to multiply images/ROIs of different sizes",3);
    }
    set_whole_cache_validity(false);
    int toff = source.mint() - this->mint();
    for (int t=this->mint(); t<=this->maxt(); t++) {
      vols[t] *= source[t + toff];
    }
    return *this;
  }


  template <class T>
  const volume4D<T>& volume4D<T>::operator/=(const volume4D<T>& source)
  {
    if (!samesize(*this,source)) {
      imthrow("Attempted to divide images/ROIs of different sizes",3);
    }
    set_whole_cache_validity(false);
    int toff = source.mint() - this->mint();
    for (int t=this->mint(); t<=this->maxt(); t++) {
      vols[t] /= source[t + toff];
    }
    return *this;
  }

  template <class T>
  volume4D<T> volume4D<T>::operator+(T num) const
  {
    set_whole_cache_validity(false);
    volume4D<T> tmp = *this;
    tmp+=num;
    return tmp;
  }

  template <class T>
  volume4D<T> volume4D<T>::operator-(T num) const
  {
    set_whole_cache_validity(false);
    volume4D<T> tmp = *this;
    tmp-=num;
    return tmp;
  }

  template <class T>
  volume4D<T> volume4D<T>::operator*(T num) const
  {
    volume4D<T> tmp = *this;
    tmp*=num;
    return tmp;
  }

  template <class T>
  volume4D<T> volume4D<T>::operator/(T num) const
  {
    set_whole_cache_validity(false);
    volume4D<T> tmp = *this;
    tmp/=num;
    return tmp;
  }

  template <class T>
  volume4D<T> volume4D<T>::operator+(const volume<T>& vol2) const
  {
    volume4D<T> tmp = *this;
    tmp+=vol2;
    return tmp;
  }

  template <class T>
  volume4D<T> volume4D<T>::operator-(const volume<T>& vol2) const
  {
    set_whole_cache_validity(false);
    volume4D<T> tmp = *this;
    tmp-=vol2;
    return tmp;
  }

  template <class T>
  volume4D<T> volume4D<T>::operator*(const volume<T>& vol2) const
  {
    volume4D<T> tmp = *this;
    tmp*=vol2;
    return tmp;
  }

  template <class T>
  volume4D<T> volume4D<T>::operator/(const volume<T>& vol2) const
  {
    set_whole_cache_validity(false);
    volume4D<T> tmp = *this;
    tmp/=vol2;
    return tmp;
  }

  template <class T>
  volume4D<T> volume4D<T>::operator+(const volume4D<T>& vol2) const
  {
    set_whole_cache_validity(false);
    volume4D<T> tmp = *this;
    tmp+=vol2;
    return tmp;
  }

  template <class T>
  volume4D<T> volume4D<T>::operator-(const volume4D<T>& vol2) const
  {
    set_whole_cache_validity(false);
    volume4D<T> tmp = *this;
    tmp-=vol2;
    return tmp;
  }

  template <class T>
  volume4D<T> volume4D<T>::operator*(const volume4D<T>& vol2) const
  {
    set_whole_cache_validity(false);
    volume4D<T> tmp = *this;
    tmp*=vol2;
    return tmp;
  }

  template <class T>
  volume4D<T> volume4D<T>::operator/(const volume4D<T>& vol2) const
  {
    set_whole_cache_validity(false);
    volume4D<T> tmp = *this;
    tmp/=vol2;
    return tmp;
  }


  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////


  // SPECIFIC INSTANTIATIONS

  // provide only these instances of the class
  //  NB: unsigned int is not included as it is not a valid AVW DT_TYPE
  template class volume<char>;
  template class volume<short>;
  template class volume<int>;
  template class volume<float>;
  template class volume<double>;
  //  template class volume<float_complex>;  // does not define less than


  template class volume4D<char>;
  template class volume4D<short>;
  template class volume4D<int>;
  template class volume4D<float>;
  template class volume4D<double>;
  //  template class volume4D<float_complex>;  // does not define less than







} // end namespace

