/*  Templated image storage class

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

#if !defined(__newimage_h)
#define __newimage_h

#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include "newmatap.h"
#include "lazy.h"
#include "lazyiterators.h"
#include "positerators.h"
#include "miscmaths/miscmaths.h"
#include "miscmaths/kernel.h"
#include "niftiio/nifti1.h"
#include "miscmaths/splinterpolator.h"
 
using namespace NEWMAT;
using namespace LAZY;
using namespace MISCMATHS;
using namespace std;

namespace NEWIMAGE {


  const bool USEASSERT=false;  // used inside imthrow to determine type of action

  void imthrow(const string& msg, int nierrnum);

  enum extrapolation { zeropad, constpad, extraslice, mirror, periodic,
                       boundsassert, boundsexception, userextrapolation };

  enum interpolation { nearestneighbour, trilinear, sinc, userkernel, 
		       userinterpolation, spline };

  enum threshtype { inclusive , exclusive };

  template <class P>
  class minmaxstuff { public: P min; P max; int minx; int miny; int minz; 
			int mint; int maxx; int maxy; int maxz; int maxt; } ;


#pragma interface  

  template <class T>
  class volume : public lazymanager {
  private:

    T* Data;
    bool data_owner;
    int SizeBound;
    int SliceOffset;

    int SlicesZ;
    int RowsY;
    int ColumnsX;

    float Xdim;
    float Ydim;
    float Zdim;

    // the following matrices map voxel coords to mm coords (nifti conventions)
    mutable Matrix StandardSpaceCoordMat;
    mutable Matrix RigidBodyCoordMat;
    mutable int StandardSpaceTypeCode;
    mutable int RigidBodyTypeCode;


    mutable int IntentCode;
    mutable float IntentParam1;
    mutable float IntentParam2;
    mutable float IntentParam3;

    mutable int SliceOrderingCode;

    mutable std::vector<int> ROIbox;
    mutable bool activeROI;
    mutable std::vector<int> Limits;

    mutable unsigned long int no_voxels;
    lazy<minmaxstuff<T>, volume<T> > minmax;
    lazy<std::vector<double>, volume<T> > sums;
    lazy<std::vector<T>, volume<T> > robustlimits;
    lazy<Matrix, volume<T> > principleaxes;
    lazy<std::vector<T>, volume<T> > percentiles;
    mutable std::vector<float> percentilepvals;
    lazy<ColumnVector, volume<T> > l_histogram;
    mutable int HISTbins;
    mutable T HISTmin;
    mutable T HISTmax;
    lazy<SPLINTERPOLATOR::Splinterpolator<T>, volume<T> > splint;

    mutable kernel interpkernel;
    mutable extrapolation p_extrapmethod;
    mutable interpolation p_interpmethod;
    mutable unsigned int splineorder;
    mutable T (*p_userextrap)(const volume<T>& , int, int, int);
    mutable float (*p_userinterp)(const volume<T>& , float, float, float);
    mutable T padvalue;
    mutable T extrapval;                 // the reference target for all extrapolations
    mutable std::vector<bool> ep_valid;  // Indicates if extrapolation can produce "valid" values.

    mutable float displayMaximum;
    mutable float displayMinimum;
    char auxFile[24];

    // Internal functions
    inline T* basicptr(int x, int y, int z)
      { return (Data + (z*RowsY + y)*ColumnsX + x); }
    inline T* basicptr(int x, int y, int z) const
      { return (Data + (z*RowsY + y)*ColumnsX + x); }
    
    int initialize(int xsize, int ysize, int zsize, T *d, bool d_owner);
    void setdefaultproperties();
    void setupsizeproperties() const;
    void setdefaultlimits() const;
    void enforcelimits(std::vector<int>& lims) const;
    void calc_no_voxels() const;
    const T& extrapolate(int x, int y, int z) const;
    float kernelinterpolation(const float x, const float y, 
			      const float z) const;
    float splineinterpolate(float x, float y, float z) const;
    float spline_interp1partial(float x, float y, float z, int dir, float *deriv) const;
    float spline_interp3partial(float x, float y, float z, float *dfdx, float *dfdy, float *dfdz) const;

    typedef T* nonsafe_fast_iterator;
    inline nonsafe_fast_iterator nsfbegin()
      { set_whole_cache_validity(false); return Data; }
    inline nonsafe_fast_iterator nsfend()
      { return Data+SizeBound; }
    template <class S, class D> friend
    void copybasicproperties(const volume<S>& source, volume<D>& dest);

    void basic_swapdimensions(int dim1, int dim2, int dim3, bool keepLRorder);
    lazy<ColumnVector, volume<T> > lazycog;  // in voxel coordinates

#ifdef EXPOSE_TREACHEROUS
  public:
#endif 
    // sampling_mat should now be avoided - use newimagevox2mm_mat instead
    bool RadiologicalFile;
    Matrix sampling_mat() const;
    void set_sform(int sform_code, const Matrix& snewmat) const;
    void set_qform(int qform_code, const Matrix& qnewmat) const;
    int  left_right_order() const; // see also newimagevox2mm_mat()
    void swapLRorder();
    void setLRorder(int LRorder);
    void makeradiological();
    void makeneurological();
    Matrix sform_mat() const { return StandardSpaceCoordMat; }
    int sform_code() const { return StandardSpaceTypeCode; }
    Matrix qform_mat() const { return RigidBodyCoordMat; }
    int qform_code() const { return RigidBodyTypeCode; }

  public:

    // CONSTRUCTORS AND DESTRUCTORS (including copy, = and reinitialize)
    volume();
    volume(const volume<T>& source);
    volume(int xsize, int ysize, int zsize);
    volume(int xsize, int ysize, int zsize, T *d, bool d_owner);
    ~volume();
    void destroy();
    const volume<T>& operator=(const volume<T>& source); 
    int reinitialize(const volume<T>& source);
    int reinitialize(int xsize, int ysize, int zsize);
    int reinitialize(int xsize, int ysize, int zsize, T *d, bool d_owner);
    int copyproperties(const volume<T>& source);
    int copydata(const volume<T>& source);

    // BASIC PROPERTIES
    inline int xsize() const { return ColumnsX; }
    inline int ysize() const { return RowsY; }
    inline int zsize() const { return SlicesZ; }
    inline float xdim() const { return Xdim; }
    inline float ydim() const { return Ydim; }
    inline float zdim() const { return Zdim; }
    inline float getDisplayMaximum() const { return displayMaximum; }
    inline float getDisplayMinimum() const { return displayMinimum; }
    inline string getAuxFile() const { return string(auxFile); }

    void setxdim(float x) { Xdim = fabs(x); }
    void setydim(float y) { Ydim = fabs(y); }
    void setzdim(float z) { Zdim = fabs(z); }
    void setdims(float x, float y, float z) 
      { setxdim(x); setydim(y); setzdim(z); }
    void setDisplayMaximumMinimum(const float maximum, const float minimum) const {  displayMaximum=maximum; displayMinimum=minimum; }
    void setDisplayMaximum(const float maximum) const { setDisplayMaximumMinimum(maximum,displayMinimum); }
    void setDisplayMinimum(const float minimum) const { setDisplayMaximumMinimum(displayMaximum,minimum); }
    void setAuxFile(const string fileName) { strncpy(auxFile,fileName.c_str(),24); }
    unsigned long int nvoxels() const { return no_voxels; }

    // ROI FUNCTIONS
    inline const std::vector<int>& limits() const { return Limits; }
    inline int limits(int n) const { return Limits[n]; }
    inline int minx() const { return Limits[0]; }
    inline int maxx() const { return Limits[3]; }
    inline int miny() const { return Limits[1]; }
    inline int maxy() const { return Limits[4]; }
    inline int minz() const { return Limits[2]; }
    inline int maxz() const { return Limits[5]; }
    inline const std::vector<int>& ROIlimits() const { return ROIbox; }
    inline int ROIlimits(int n) const { return ROIbox[n]; }
    inline bool usingROI() const { return activeROI; }
    void setROIlimits(int x0, int y0, int z0, int x1, int y1, int z1) const;
    void setROIlimits(const std::vector<int>& lims) const;
    void activateROI() const; 
    void deactivateROI() const;
    volume<T> ROI() const;  // returns a new volume = ROI
    int copyROIonly(const volume<T>& source);

    // Volume<->ColumnVector conversion

    ReturnMatrix vec(const volume<T>& mask) const;
    ReturnMatrix vec() const;
    void insert_vec(const ColumnVector& pvec, const volume<T>& pmask);
    void insert_vec(const ColumnVector& pvec);
    vector<int> labelToCoord(const long label) const;

    // SECONDARY PROPERTIES
    // maps *NEWIMAGE* voxel coordinates to mm (consistent with FSLView mm)
    // NB: do not try to determine left-right order from this matrix
    // sampling_mat should now be avoided - use newimagevox2mm_mat instead
    Matrix newimagevox2mm_mat() const;
    Matrix niftivox2newimagevox_mat() const;

    int intent_code() const { return IntentCode; }
    float intent_param(int n) const;
    void set_intent(int intent_code, float p1, float p2, float p3) const;

    T min() const { return minmax().min; }
    T max() const { return minmax().max; }
    int mincoordx() const { return minmax().minx; }
    int mincoordy() const { return minmax().miny; }
    int mincoordz() const { return minmax().minz; }
    int maxcoordx() const { return minmax().maxx; }
    int maxcoordy() const { return minmax().maxy; }
    int maxcoordz() const { return minmax().maxz; }
    double sum() const { return sums()[0]; }
    double sumsquares() const { return sums()[1]; }
    double mean() const { return sum()/((double) no_voxels); }
    double variance() const { double n=(double) no_voxels; 
		return (n/(n-1))*(sumsquares()/n - mean()*mean()); }
    double stddev() const { return sqrt(variance()); }
    T robustmin() const { return robustlimits()[0]; }
    T robustmax() const { return robustlimits()[1]; }
    ColumnVector principleaxis(int n) const;
    Matrix principleaxes_mat() const;
    T percentile(float pvalue) const;  // argument in range [0.0 , 1.0]
    std::vector<float> percentilepvalues() const { return percentilepvals; }
    ColumnVector histogram(int nbins) const;
    ColumnVector histogram(int nbins, T minval, T maxval) const;
    int histbins() const { return HISTbins; }
    T histmin() const { return HISTmin; }
    T histmax() const { return HISTmax; }

    lazy<T, volume<T> > backgroundval;
    ColumnVector cog(const string& coordtype="voxel") const;

    // SECONDARY PROPERTIES (using mask)
    T min(const volume<T>& mask) const;
    T max(const volume<T>& mask) const;
    int mincoordx(const volume<T>& mask) const;
    int mincoordy(const volume<T>& mask) const;
    int mincoordz(const volume<T>& mask) const;
    int maxcoordx(const volume<T>& mask) const;
    int maxcoordy(const volume<T>& mask) const;
    int maxcoordz(const volume<T>& mask) const;
    double sum(const volume<T>& mask) const;
    double sumsquares(const volume<T>& mask) const;
    double mean(const volume<T>& mask) const;
    double variance(const volume<T>& mask) const;
    double stddev(const volume<T>& mask) const { return sqrt(variance(mask)); }
    T robustmin(const volume<T>& mask) const;
    T robustmax(const volume<T>& mask) const;
    T percentile(float pvalue, const volume<T>& mask) const;  // arg in [0,1]
    ColumnVector histogram(int nbins, const volume<T>& mask) const;
    ColumnVector histogram(int nbins, T minval, T maxval, const volume<T>& mask) 
      const;


    // DATA ACCESS FUNCTIONS (iterators)
    typedef const T* fast_const_iterator;

    inline fast_const_iterator fbegin() const 
      { return Data; }
    inline fast_const_iterator fend() const 
      { return Data+SizeBound; }


    // BASIC DATA ACCESS FUNCTIONS
    inline bool in_bounds(int x, int y, int z) const
      { return ( (x>=0) && (y>=0) && (z>=0) 
		 && (x<ColumnsX) && (y<RowsY) && (z<SlicesZ) ); }
    bool in_bounds(float x, float y, float z) const
    {
      int ix=((int) floor(x)); 
      int iy=((int) floor(y)); 
      int iz=((int) floor(z));
      return((ix>=0) && (iy>=0) && (iz>=0) && ((ix+1)<ColumnsX) && ((iy+1)<RowsY) && ((iz+1)<SlicesZ));
    }
    bool in_extraslice_bounds(float x, float y, float z) const
    {
      int ix=((int) floor(x)); 
      int iy=((int) floor(y)); 
      int iz=((int) floor(z));
      return((ix>=-1) && (iy>=-1) && (iz>=-1) && (ix<ColumnsX) && (iy<RowsY) && (iz<SlicesZ));
    }
    inline bool valid(int x, int y, int z) const 
    { 
      return((ep_valid[0] || (x>=0 && x<ColumnsX)) && (ep_valid[1] || (y>=0 && y<RowsY)) && (ep_valid[2] || (z>=0 && z<SlicesZ)));
    }
    bool valid(float x, float y, float z) const
    {
      int ix=((int) floor(x)); 
      int iy=((int) floor(y)); 
      int iz=((int) floor(z));
      return((ep_valid[0] || (ix>=0 && (ix+1)<ColumnsX)) && (ep_valid[1] || (iy>=0 && (iy+1)<RowsY)) && (ep_valid[2] || (iz>=0 && (iz+1)<SlicesZ)));
    }
    inline T& operator()(int x, int y, int z)
      { set_whole_cache_validity(false); 
        if (in_bounds(x,y,z)) return *(basicptr(x,y,z)); 
	else                  return const_cast<T& > (extrapolate(x,y,z)); }
    inline const T& operator()(int x, int y, int z)
      const { 	if (in_bounds(x,y,z)) return *(basicptr(x,y,z)); 
	        else                  return extrapolate(x,y,z); }
    float interpolate(float x, float y, float z) const;
    float interpolate(float x, float y, float z, bool *ep) const;
    float interp1partial(// Input
                         float x, float y, float z,  // Co-ordinates to get value for
                         int     dir,                   // Direction for partial, 0->x, 1->y, 2->z
                         // Output
                         float  *pderiv                // Derivative returned here
                         ) const;
    float interp3partial(// Input
                         float x, float y, float z,             // Co-ordinate to get value for
                         // Output
                         float *dfdx, float *dfdy, float *dfdz  // Partials
                         ) const;

    inline T& value(int x, int y, int z)
      { set_whole_cache_validity(false); 
        return *(basicptr(x,y,z)); }
    inline const T& value(int x, int y, int z) const
      { return *(basicptr(x,y,z)); }
    float interpolatevalue(float x, float y, float z) const;


    // SECONDARY FUNCTIONS
    void setextrapolationmethod(extrapolation extrapmethod) const { p_extrapmethod = extrapmethod; }
    extrapolation getextrapolationmethod() const { return(p_extrapmethod); }
    void setpadvalue(T padval) const { padvalue = padval; }
    T getpadvalue() const { return padvalue; }
    void defineuserextrapolation(T (*extrap)(
             const volume<T>& , int, int, int)) const;

    void setinterpolationmethod(interpolation interpmethod) const;
    interpolation getinterpolationmethod() const { return p_interpmethod; }
    void setsplineorder(unsigned int order) const;
    unsigned int getsplineorder() const { return(splineorder); }
    void setextrapolationvalidity(bool xv, bool yv, bool zv) const { ep_valid[0]=xv; ep_valid[1]=yv; ep_valid[2]=zv; }
    std::vector<bool> getextrapolationvalidity() const { return(ep_valid); }
    void defineuserinterpolation(float (*interp)(
             const volume<T>& , float, float, float)) const;
    void definekernelinterpolation(const ColumnVector& kx, 
				   const ColumnVector& ky,
				   const ColumnVector& kz, 
				   int wx, int wy, int wz) const;  // full-width
    void definekernelinterpolation(const volume<T>& vol) const;
    void definesincinterpolation(const string& sincwindowtype,
				 int w, int nstore=1201) const;  // full-width
    void definesincinterpolation(const string& sincwindowtype,
				 int wx, int wy, int wz, int nstore=1201) const;
                                  // full-width

    inline void getneighbours(int x, int y, int z, 
			      T &v000, T &v001, T &v010,
			      T &v011, T &v100, T &v101,
			      T &v110, T &v111) const;
    inline void getneighbours(int x, int y, int z, 
			      T &v000, T &v010,
			      T &v100, T &v110) const;
    
    

    // ARITHMETIC FUNCTIONS
    T operator=(T val); 
    const volume<T>& operator+=(T val); 
    const volume<T>& operator-=(T val); 
    const volume<T>& operator*=(T val); 
    const volume<T>& operator/=(T val); 
    const volume<T>& operator+=(const volume<T>& source); 
    const volume<T>& operator-=(const volume<T>& source); 
    const volume<T>& operator*=(const volume<T>& source); 
    const volume<T>& operator/=(const volume<T>& source); 

    volume<T> operator+(T num) const;
    volume<T> operator-(T num) const;
    volume<T> operator*(T num) const;
    volume<T> operator/(T num) const;
    volume<T> operator+(const volume<T>& vol2) const;
    volume<T> operator-(const volume<T>& vol2) const;
    volume<T> operator*(const volume<T>& vol2) const;
    volume<T> operator/(const volume<T>& vol2) const;

    template <class S>
    friend volume<S> operator+(S num, const volume<S>& vol);
    template <class S>
    friend volume<S> operator-(S num, const volume<S>& vol);
    template <class S>
    friend volume<S> operator*(S num, const volume<S>& vol);
    template <class S>
    friend volume<S> operator/(S num, const volume<S>& vol);
    template <class S>
    friend volume<S> operator-(const volume<S>& vol);

    // Comparisons. These are used for "spatial" purposes
    // so that if data is identical and all the "spatial
    // fields of the header are identical then the volumes
    // are considered identical.
    template <class S>
    friend bool operator==(const volume<S>& v1, const volume<S>& v2);
    template <class S>
    friend bool operator!=(const volume<S>& v1, const volume<S>& v2); // { return(!(v1==v2)); }
    
    // GENERAL MANIPULATION

    void binarise(T lowerth, T upperth, threshtype tt=inclusive);
    void binarise(T thresh) { this->binarise(thresh,this->max(),inclusive); }
    void threshold(T lowerth, T upperth, threshtype tt=inclusive);
    void threshold(T thresh) { this->threshold(thresh,this->max(),inclusive); }
    // valid entries for dims are +/- 1, 2, 3 (and for newx, etc they are x, -x, y, -y, z, -z)
    void swapdimensions(int dim1, int dim2, int dim3);
    void swapdimensions(const string& newx, const string& newy, const string& newz);
    Matrix swapmat(int dim1, int dim2, int dim3) const;
    Matrix swapmat(const string& newx, const string& newy, const string& newz) const;

     
    // CONVERSION FUNCTIONS
    template <class S, class D> friend
    void copyconvert(const volume<S>& source, volume<D>& dest);
      
  };


  // HELPER FUNCTIONS

  template <class T>
  long int no_mask_voxels(const volume<T>& mask);

  template <class S, class D>
  void convertbuffer(const S* source, D* dest, int len);

  template <class S, class D>
  void convertbuffer(const S* source, D* dest, int len, float slope, float intercept);

  template <class S1, class S2>
  bool samesize(const volume<S1>& vol1, const volume<S2>& vol2, bool checkdim=false);

  template <class S1, class S2>
  bool sameabssize(const volume<S1>& vol1, const volume<S2>& vol2, bool checkdim=false);

  template <class S1, class S2>
  bool samedim(const volume<S1>& vol1, const volume<S2>& vol2);


  //////////////////////// VOLUME4D CLASS /////////////////////////////


  template <class T>
  class volume4D : public lazymanager {
  private:
    std::vector<volume<T> > vols;
    float p_TR;

    mutable std::vector<int> ROIbox;
    mutable bool activeROI;
    mutable std::vector<int> Limits;

    mutable extrapolation p_extrapmethod;
    mutable interpolation p_interpmethod;
    mutable T (*p_userextrap)(const volume<T>& , int, int, int);
    mutable float (*p_userinterp)(const volume<T>& , float, float, float);
    mutable T p_padval;

    lazy<minmaxstuff<T>, volume4D<T> > tsminmax;
    lazy<std::vector<double>, volume4D<T> > sums;
    lazy<std::vector<T>, volume4D<T> > robustlimits;
    lazy<std::vector<T>, volume4D<T> > percentiles;
    mutable std::vector<float> percentilepvals;
    lazy<ColumnVector, volume4D<T> > l_histogram;
    mutable int HISTbins;
    mutable T HISTmin;
    mutable T HISTmax;

    // Internal functions
    int initialize(int xsize, int ysize, int zsize, int tsize, T *d=0);
    void setdefaultproperties();
    void setdefaultlimits() const;
    void enforcelimits(std::vector<int>& lims) const;
    template <class S, class D> friend
    void copybasicproperties(const volume<S>& source, volume4D<D>& dest);
    template <class S, class D> friend
    void copybasicproperties(const volume4D<S>& source, volume<D>& dest);
    template <class S, class D> friend
    void copybasicproperties(const volume4D<S>& source, volume4D<D>& dest);

#ifdef EXPOSE_TREACHEROUS
  public:
#endif 
    // sampling_mat should now be avoided - use newimagevox2mm_mat instead
    Matrix sampling_mat() const;
    void set_sform(int sform_code, const Matrix& snewmat) const;
    void set_qform(int qform_code, const Matrix& qnewmat) const;
    int  left_right_order() const;
    void swapLRorder();
    void setLRorder(int LRorder);
    void makeradiological();
    void makeneurological();
    Matrix sform_mat() const;
    int sform_code() const;
    Matrix qform_mat() const;
    int qform_code() const;

  public:
    // CONSTRUCTORS AND DESTRUCTORS (including copy, = and reinitialize)
    volume4D();
    volume4D(const volume4D<T>& source);
    volume4D(int xsize, int ysize, int zsize, int tsize, T *d=0);
    ~volume4D();
    void destroy();
    const volume4D<T>& operator=(const volume4D<T>& source); 
    const volume4D<T>& operator=(const volume<T>& source); 
    int reinitialize(const volume4D<T>& source);
    int reinitialize(int xsize, int ysize, int zsize, int tsize, T *d=0);
    int copyproperties(const volume4D<T>& source);
    int copyproperties(const volume<T>& source);
    int copyvolumes(const volume4D<T>& source);
    
    // DATA ACCESS
    inline bool in_bounds(int x, int y, int z) const
      { return ( (vols.size()>0) && vols[0].in_bounds(x,y,z) ); }
    bool in_bounds(float x, float y, float z) const
      { return ( (vols.size()>0) && vols[0].in_bounds(x,y,z) ); }
    inline bool in_bounds(int t) const
      { return ( (t>=0) && (t<this->tsize()) ); }
    inline bool in_bounds(int x, int y, int z, int t) const
      { return ( (t>=0) && (t<this->tsize()) && 
		 vols[Limits[3]].in_bounds(x,y,z) ); }
    bool in_bounds(float x, float y, float z, int t) const
      { return ( (t>=0) && (t<this->tsize()) && 
		 vols[Limits[3]].in_bounds(x,y,z) ); }
    bool valid(int x, int y, int z) const {return(vols.size()>0 && vols[0].valid(x,y,z));}
    bool valid(float x, float y, float z) const {return(vols.size()>0 && vols[0].valid(x,y,z));}

    inline T& operator()(int x, int y, int z, int t) 
      { set_whole_cache_validity(false);  
        if (!in_bounds(t)) imthrow("Out of Bounds (time index)",5); 
	return vols[t](x,y,z); }
    inline const T& operator()(int x, int y, int z, int t) const 
      { if (!in_bounds(t)) imthrow("Out of Bounds (time index)",5); 
        return vols[t](x,y,z); }

    inline T& value(int x, int y, int z, int t) 
      { set_whole_cache_validity(false);  return vols[t].value(x,y,z); }
    inline const T& value(int x, int y, int z, int t) const 
      { return vols[t].value(x,y,z); }

    inline volume<T>& operator[](int t) {
      set_whole_cache_validity(false);  
      if (!in_bounds(t)) imthrow("Out of Bounds (time index)",5);   
      return vols[t]; }
    inline const volume<T>& operator[](int t) const 
      { if (!in_bounds(t)) imthrow("Out of Bounds (time index)",5);  
        return vols[t]; }
    
    // SINGLE VOXEL TIME-SERIES ACCESS
    ReturnMatrix voxelts(int x, int y, int z) const;
    void setvoxelts(const ColumnVector& ts, int x, int y, int z);
    
    // WHOLE VOLUME MODIFIERS
    void addvolume(const volume<T>& source);
    void addvolume(const volume4D<T>& source);
    void insertvolume(const volume<T>& source, int t);
    void deletevolume(int t);
    void clear();  // deletes all volumes

    // MATRIX <-> VOLUME4D CONVERSIONS
    ReturnMatrix matrix(const volume<T>& mask) const;
    ReturnMatrix matrix(const volume<T>& mask, vector<long>& voxelLabels) const;
    ReturnMatrix matrix() const;
    void setmatrix(const Matrix& newmatrix, const volume<T>& mask, 
		   const T pad=0);
    void setmatrix(const Matrix& newmatrix); 
    volume<int> vol2matrixkey(volume<T>& mask); //returns a volume with numbers in relating to matrix colnumbers
    ReturnMatrix matrix2volkey(volume<T>& mask);


    // SIZE AND DIMENSIONS
    inline int xsize() const 
      { if (vols.size()>0) return vols[0].xsize(); else return 0; }
    inline int ysize() const
      { if (vols.size()>0) return vols[0].ysize(); else return 0; }
    inline int zsize() const
      { if (vols.size()>0) return vols[0].zsize(); else return 0; }
    inline int tsize() const { return (int) vols.size(); }
    inline float xdim() const
      { if (vols.size()>0) return vols[0].xdim(); else return 1.0; }
    inline float ydim() const
      { if (vols.size()>0) return vols[0].ydim(); else return 1.0; }
    inline float zdim() const
      { if (vols.size()>0) return vols[0].zdim(); else return 1.0; }
    inline float tdim() const { return p_TR; }
    inline float TR() const { return p_TR; }
    inline float getDisplayMaximum() const { if (vols.size()>0) return vols[0].getDisplayMaximum(); else return 0; }
    inline float getDisplayMinimum() const { if (vols.size()>0) return vols[0].getDisplayMinimum(); else return 0; }
    inline string getAuxFile() const { if (vols.size()>0) return vols[0].getAuxFile(); else return string(""); }

    void setxdim(float x);
    void setydim(float y);
    void setzdim(float z);
    void settdim(float tr) { p_TR = fabs(tr); }
    void setTR(float tr)   { settdim(tr); }
    void setdims(float x, float y, float z, float tr) 
      { setxdim(x); setydim(y); setzdim(z); settdim(tr); }
    unsigned long int nvoxels() const 
      { if (vols.size()>0) return vols[0].nvoxels(); else return 0; }
    int ntimepoints() const { return this->tsize(); }
    void setDisplayMaximumMinimum(const float maximum, const float minimum) const;
    void setDisplayMaximum(const float maximum) const { setDisplayMaximumMinimum(maximum,vols[0].getDisplayMinimum()); }
    void setDisplayMinimum(const float minimum) const { setDisplayMaximumMinimum(vols[0].getDisplayMaximum(),minimum); }
    void setAuxFile(const string fileName);

    void setinterpolationmethod(interpolation interpmethod) const;
    interpolation getinterpolationmethod() const;
    void setsplineorder(unsigned int order) const;
    unsigned int getsplineorder() const;
    void setextrapolationvalidity(bool xv, bool yv, bool zv) const;
    std::vector<bool> getextrapolationvalidity() const;
    void setextrapolationmethod(extrapolation extrapmethod) const;
    extrapolation getextrapolationmethod() const;
    void setpadvalue(T padval) const;
    T getpadvalue() const;
    void defineuserinterpolation(float (*interp)(
		       const volume<T>& , float, float, float)) const;
    void defineuserextrapolation(T (*extrap)(
		       const volume<T>& , int, int, int)) const;
    void definekernelinterpolation(const ColumnVector& kx, 
				   const ColumnVector& ky,
				   const ColumnVector& kz, 
				   int wx, int wy, int wz) const;
    void definekernelinterpolation(const volume4D<T>& vol) const;
    void definekernelinterpolation(const volume<T>& vol) const;
    void definesincinterpolation(const string& sincwindowtype,
				 int w, int nstore=1201) const;
    void definesincinterpolation(const string& sincwindowtype,
				 int wx, int wy, int wz, int nstore=1201) const;

    // ROI FUNCTIONS
    inline const std::vector<int>& limits() const { return Limits; }
    inline int limits(int n) const { return Limits[n]; }
    inline int minx() const { return Limits[0]; }
    inline int maxx() const { return Limits[4]; }
    inline int miny() const { return Limits[1]; }
    inline int maxy() const { return Limits[5]; }
    inline int minz() const { return Limits[2]; }
    inline int maxz() const { return Limits[6]; }
    inline int mint() const { return Limits[3]; }
    inline int maxt() const { return Limits[7]; }
    inline const std::vector<int>& ROIlimits() const { return ROIbox; }
    inline int ROIlimits(int n) const { return ROIbox[n]; }
    inline bool usingROI() const { return activeROI; }
    void setROIlimits(int x0, int y0, int z0, int t0,
                      int x1, int y1, int z1, int t1) const;
    void setROIlimits(int t0, int t1) const;
    void setROIlimits(int x0, int y0, int z0, int x1, int y1, int z1) const; 
    void setROIlimits(const std::vector<int>& lims) const;
    void activateROI() const; 
    void deactivateROI() const;
    volume4D<T> ROI() const;  // returns a new volume = ROI
    int copyROIonly(const volume4D<T>& source);


    // SECONDARY PROPERTIES
    // maps *NEWIMAGE* voxel coordinates to mm (consistent with FSLView mm)
    // NB: do not try to determine left-right order from this matrix
    // sampling_mat should now be avoided - use newimagevox2mm_mat instead
    Matrix newimagevox2mm_mat() const;
    Matrix niftivox2newimagevox_mat() const;
    
    int intent_code() const;
    float intent_param(int n) const;
    void set_intent(int intent_code, float p1, float p2, float p3) const;

    T min() const { return tsminmax().min; }
    T max() const { return tsminmax().max; }
    int mincoordx() const { return tsminmax().minx; }
    int mincoordy() const { return tsminmax().miny; }
    int mincoordz() const { return tsminmax().minz; }
    int mincoordt() const { return tsminmax().mint; }
    int maxcoordx() const { return tsminmax().maxx; }
    int maxcoordy() const { return tsminmax().maxy; }
    int maxcoordz() const { return tsminmax().maxz; }
    int maxcoordt() const { return tsminmax().maxt; }
    double sum() const { return sums()[0]; }
    double sumsquares() const { return sums()[1]; }
    double mean() const { return sum()/(Max(1.0,(double) nvoxels()*ntimepoints()));}
    double variance() const { double n=(double) nvoxels() * ntimepoints(); 
		return (n/(n-1))*(sumsquares()/n - mean()*mean()); }
    double stddev() const { return sqrt(variance()); }
    T robustmin() const { return robustlimits()[0]; }
    T robustmax() const { return robustlimits()[1]; }
    Matrix principleaxes_mat() const;
    T percentile(float pvalue) const;  // argument in range [0.0 , 1.0]
    std::vector<float> percentilepvalues() const { return percentilepvals; }
    ColumnVector histogram(int nbins) const;
    ColumnVector histogram(int nbins, T minval, T maxval) const;
    int histbins() const { return HISTbins; }
    T histmin() const { return HISTmin; }
    T histmax() const { return HISTmax; }

    // SECONDARY PROPERTIES (using a 3D mask)
    T min(const volume<T>& mask) const;
    T max(const volume<T>& mask) const;
    int mincoordx(const volume<T>& mask) const;
    int mincoordy(const volume<T>& mask) const;
    int mincoordz(const volume<T>& mask) const;
    int maxcoordx(const volume<T>& mask) const;
    int maxcoordy(const volume<T>& mask) const;
    int maxcoordz(const volume<T>& mask) const;
    double sum(const volume<T>& mask) const;
    double sumsquares(const volume<T>& mask) const;
    double mean(const volume<T>& mask) const;
    double variance(const volume<T>& mask) const;
    double stddev(const volume<T>& mask) const { return sqrt(variance(mask)); }
    T percentile(float pvalue, const volume<T>& mask) const;  // arg in [0,1]
    T robustmin(const volume<T>& mask) const;
    T robustmax(const volume<T>& mask) const;
    ColumnVector histogram(int nbins, const volume<T>& mask) const;
    ColumnVector histogram(int nbins, T minval, T maxval, const volume<T>& mask)
      const;

    // SECONDARY PROPERTIES (using a 4D mask)
    T min(const volume4D<T>& mask) const;
    T max(const volume4D<T>& mask) const;
    int mincoordx(const volume4D<T>& mask) const;
    int mincoordy(const volume4D<T>& mask) const;
    int mincoordz(const volume4D<T>& mask) const;
    int maxcoordx(const volume4D<T>& mask) const;
    int maxcoordy(const volume4D<T>& mask) const;
    int maxcoordz(const volume4D<T>& mask) const;
    double sum(const volume4D<T>& mask) const;
    double sumsquares(const volume4D<T>& mask) const;
    double mean(const volume4D<T>& mask) const;
    double variance(const volume4D<T>& mask) const;
    double stddev(const volume4D<T>& mask) const { return sqrt(variance(mask));}
    T percentile(float pvalue, const volume4D<T>& mask) const;  // arg in [0,1]
    T robustmin(const volume4D<T>& mask) const;
    T robustmax(const volume4D<T>& mask) const;
    ColumnVector histogram(int nbins, const volume4D<T>& mask) const;
    ColumnVector histogram(int nbins, T minval, T maxval, const volume4D<T>& mask)
      const;


    // GENERAL MANIPULATION
    void binarise(T lowerth, T upperth, threshtype tt=inclusive);
    void binarise(T thresh) { this->binarise(thresh,this->max(),inclusive); }
    void threshold(T lowerth, T upperth, threshtype tt=inclusive);
    void threshold(T thresh) { this->threshold(thresh,this->max(),inclusive); }
    // valid entries for dims are +/- 1, 2, 3 (and for newx, etc they are x, -x, y, -y, z, -z)
    void swapdimensions(int dim1, int dim2, int dim3);
    void swapdimensions(const string& newx, const string& newy, const string& newz);
    Matrix swapmat(int dim1, int dim2, int dim3) const;
    Matrix swapmat(const string& newx, const string& newy, const string& newz) const;


    // ARITHMETIC FUNCTIONS
    T operator=(T val); 
    const volume4D<T>& operator+=(T val); 
    const volume4D<T>& operator-=(T val); 
    const volume4D<T>& operator*=(T val); 
    const volume4D<T>& operator/=(T val); 
    const volume4D<T>& operator+=(const volume<T>& source); 
    const volume4D<T>& operator-=(const volume<T>& source); 
    const volume4D<T>& operator*=(const volume<T>& source); 
    const volume4D<T>& operator/=(const volume<T>& source); 
    const volume4D<T>& operator+=(const volume4D<T>& source); 
    const volume4D<T>& operator-=(const volume4D<T>& source); 
    const volume4D<T>& operator*=(const volume4D<T>& source); 
    const volume4D<T>& operator/=(const volume4D<T>& source); 

    volume4D<T> operator+(T num) const;
    volume4D<T> operator-(T num) const;
    volume4D<T> operator*(T num) const;
    volume4D<T> operator/(T num) const;
    volume4D<T> operator+(const volume<T>& vol2) const;
    volume4D<T> operator-(const volume<T>& vol2) const;
    volume4D<T> operator*(const volume<T>& vol2) const;
    volume4D<T> operator/(const volume<T>& vol2) const;
    volume4D<T> operator+(const volume4D<T>& vol2) const;
    volume4D<T> operator-(const volume4D<T>& vol2) const;
    volume4D<T> operator*(const volume4D<T>& vol2) const;
    volume4D<T> operator/(const volume4D<T>& vol2) const;

    template <class S>
    friend volume4D<S> operator+(S num, const volume4D<S>& vol);
    template <class S>
    friend volume4D<S> operator-(S num, const volume4D<S>& vol);
    template <class S>
    friend volume4D<S> operator*(S num, const volume4D<S>& vol);
    template <class S>
    friend volume4D<S> operator/(S num, const volume4D<S>& vol);
    template <class S>
    friend volume4D<S> operator-(const volume4D<S>& vol);
    template <class S>
    friend volume4D<S> operator+(const volume<S>& v1, const volume4D<S>& v2);
    template <class S>
    friend volume4D<S> operator-(const volume<S>& v1, const volume4D<S>& v2);
    template <class S>
    friend volume4D<S> operator*(const volume<S>& v1, const volume4D<S>& v2);
     
    // CONVERSION FUNCTIONS
    template <class S, class D> friend
    void copyconvert(const volume4D<S>& source, volume4D<D>& dest);
  };


  // HELPER FUNCTIONS

  template <class S1, class S2>
  bool samesize(const volume4D<S1>& vol1, const volume4D<S2>& vol2, bool checkdim=false);

  template <class S1, class S2>
  bool sameabssize(const volume4D<S1>& vol1, const volume4D<S2>& vol2, bool checkdim=false);

  template <class S1, class S2>
  bool samedim(const volume4D<S1>& vol1, const volume4D<S2>& vol2);

////////////////////////////////////////////////////////////////////////
///////////////////////// DEBUGGING FUNCTIONS //////////////////////////
////////////////////////////////////////////////////////////////////////

  template <class T>
  void print_info(const volume<T>& source, const string& name)
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
  void print_info(const volume4D<T>& source, const string& name)
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


////////////////////////////////////////////////////////////////////////
///////////////////////// INLINE DEFINITIONS ///////////////////////////
////////////////////////////////////////////////////////////////////////

  template <class T>
  inline void volume<T>::getneighbours(int x, int y, int z, 
					    T &v000, T &v001, T &v010,
					    T &v011, T &v100, T &v101,
					    T &v110, T &v111) const {
    T *ptr = basicptr(x,y,z);
    v000 = *ptr;
    ptr++; 
    v100 = *ptr;
    ptr+= ColumnsX;  
    v110 = *ptr;
    ptr--; 
    v010 = *ptr;
    ptr += SliceOffset;
    v011 = *ptr;
    ptr++;
    v111 = *ptr;
    ptr-= ColumnsX;
    v101 = *ptr;
    ptr--;
    v001 = *ptr;
  }

	
  template <class T>
  inline void volume<T>::getneighbours(int x, int y, int z, 
					    T &v000, T &v010, 
					    T &v100, T &v110) const {
    T *ptr = basicptr(x,y,z);
    v000 = *ptr;
    ptr++; 
    v100 = *ptr;
    ptr+= ColumnsX;  
    v110 = *ptr;
    ptr--; 
    v010 = *ptr;
  }
	

//    template <class T>
//    inline volume<T>::iterator volume<T>::begin()
//    { return volume<T>::iterator(basicptr(Limits[0],Limits[1],Limits[2]),
//  		    (lazymanager*) this, 
//  		    Limits[0],Limits[1],Limits[2],
//  		    Limits[0],Limits[1],Limits[2],
//  		    Limits[3],Limits[4],Limits[5],
//  		    ColumnsX,SliceOffset); }

//    template <class T>
//    inline volume<T>::iterator volume<T>::end()
//    { volume<T>::iterator tmp(basicptr(Limits[3],Limits[4],Limits[5]),
//  		 (lazymanager*) this,
//  		 Limits[3],Limits[4],Limits[5],
//  		 Limits[0],Limits[1],Limits[2],
//  		 Limits[3],Limits[4],Limits[5],
//  		 ColumnsX,SliceOffset); 
//                   return ++tmp; }

//    template <class T>
//    inline volume<T>::const_iterator volume<T>::begin() const 
//    { return 
//        volume<T>::const_iterator(basicptr(Limits[0],Limits[1],Limits[2]),
//  			  Limits[0],Limits[1],Limits[2],
//  			  Limits[0],Limits[1],Limits[2],
//  			  Limits[3],Limits[4],Limits[5],
//  			  ColumnsX,SliceOffset); }

//    template <class T>
//    inline volume<T>::const_iterator volume<T>::end() const 
//    { volume<T>::const_iterator tmp(basicptr(Limits[3],Limits[4],Limits[5]),
//  		       Limits[3],Limits[4],Limits[5],
//  		       Limits[0],Limits[1],Limits[2],
//  		       Limits[3],Limits[4],Limits[5],
//  		       ColumnsX,SliceOffset); 
//                         return ++tmp; }
  

	     
////////////////////////////////////////////////////////////////////////
/////////////////////////// HELPER FUNCTIONS ///////////////////////////
////////////////////////////////////////////////////////////////////////


  template <class T>
  long int no_mask_voxels(const volume<T>& mask) 
  {
     long int n=0;
     for (int z=mask.minz(); z<=mask.maxz(); z++) {
       for (int y=mask.miny(); y<=mask.maxy(); y++) {
         for (int x=mask.minx(); x<=mask.maxx(); x++) {
	   if (mask.value(x,y,z)>(T) 0.5) n++;
         }
       }
     }
     return n;
  }
  
  template <class T>
  long int no_mask_voxels(const volume4D<T>& mask) 
  {
     long int n=0;
     for (int t=mask.mint(); t<=mask.maxt(); t++) {
       for (int z=mask.minz(); z<=mask.maxz(); z++) {
	 for (int y=mask.miny(); y<=mask.maxy(); y++) {
	   for (int x=mask.minx(); x<=mask.maxx(); x++) {
	     if (mask.value(x,y,z,t)>(T) 0.5) n++;
	   }
	 }
       }
     }
     return n;
  }

  template <class T>
  long int no_timepoints(const volume4D<T>& vol) 
  {
     return (vol.maxt()-vol.mint()+1);
  }

  template <class S, class D>
  void convertbuffer(const S* source, D* dest, int len)
  {
    D* dptr=dest;
    for (const S* sptr=source; sptr<(source+len); sptr++) {
      *dptr = (D) *sptr;
      dptr++;
    }
  }

  template <class S, class D>
  void convertbuffer(const S* source, D* dest, int len, float slope, float intercept)
  {
    D* dptr=dest;
    for (const S* sptr=source; sptr<(source+len); sptr++) {
      *dptr = (D) ((*sptr) * slope + intercept);
      dptr++;
    }
  }


  template <class S1, class S2>
  bool samesize(const volume<S1>& vol1, const volume<S2>& vol2, bool checkdim)
  {  
    bool same = ( ( (vol1.maxx()-vol1.minx())==(vol2.maxx()-vol2.minx()) ) && 
	          ( (vol1.maxy()-vol1.miny())==(vol2.maxy()-vol2.miny()) ) && 
	          ( (vol1.maxz()-vol1.minz())==(vol2.maxz()-vol2.minz()) ) );
    if (checkdim) same = same && samedim(vol1,vol2);
    return(same); 
  }
  
  template <class S1, class S2>
  bool samesize(const volume4D<S1>& vol1, const volume4D<S2>& vol2, bool checkdim)
  {  
    bool same = ( ( (vol1.maxt()-vol1.mint())==(vol2.maxt()-vol2.mint()) ) && 
	          ( (vol1.tsize()<=0) || (vol2.tsize()<=0) || samesize(vol1[0],vol2[0]) ) );
    if (checkdim) same = same && samedim(vol1,vol2);
    return(same);
  }
  
  template <class S1, class S2>
  bool sameabssize(const volume<S1>& vol1, const volume<S2>& vol2, bool checkdim)
  {  
    bool same = ( ( (vol1.xsize())==(vol2.xsize()) ) && 
	          ( (vol1.ysize())==(vol2.ysize()) ) && 
	          ( (vol1.zsize())==(vol2.zsize()) ) );
    if (checkdim) same = same && samedim(vol1,vol2);
    return(same); 
  }
  
  template <class S1, class S2>
  bool sameabssize(const volume4D<S1>& vol1, const volume4D<S2>& vol2, bool checkdim)
  {  
    bool same = ( ( (vol1.tsize())==(vol2.tsize()) ) && 
	          ( (vol1.tsize()==0) || samesize(vol1[0],vol2[0]) ) );
    if (checkdim) same = same && samedim(vol1,vol2);
    return(same);
  }
  
  template <class S1, class S2>
  bool samedim(const volume<S1>& vol1, const volume<S2>& vol2)
  {
    return(std::abs(vol1.xdim()-vol2.xdim())<1e-6 && 
           std::abs(vol1.ydim()-vol2.ydim())<1e-6 && 
           std::abs(vol1.zdim()-vol2.zdim())<1e-6);
  }

  template <class S1, class S2>
  bool samedim(const volume4D<S1>& vol1, const volume4D<S2>& vol2)
  {
    return(std::abs(vol1.tdim()-vol2.tdim())<1e-6 && samedim(vol1[0],vol2[0]));
  }

////////////////////////////////////////////////////////////////////////
/////////////////////////// FRIEND FUNCTIONS ///////////////////////////
////////////////////////////////////////////////////////////////////////

  template <class S, class D>
  void copybasicproperties(const volume<S>& source, volume<D>& dest)
  {
    // set up properties (except lazy ones)
    dest.Xdim = source.Xdim;
    dest.Ydim = source.Ydim;
    dest.Zdim = source.Zdim;
    
    dest.StandardSpaceCoordMat = source.StandardSpaceCoordMat;
    dest.RigidBodyCoordMat = source.RigidBodyCoordMat;
    dest.StandardSpaceTypeCode = source.StandardSpaceTypeCode;
    dest.RigidBodyTypeCode = source.RigidBodyTypeCode;

    dest.RadiologicalFile = source.RadiologicalFile;

    dest.IntentCode = source.IntentCode;
    dest.IntentParam1 = source.IntentParam1;
    dest.IntentParam2 = source.IntentParam2;
    dest.IntentParam3 = source.IntentParam3;

    dest.SliceOrderingCode = source.SliceOrderingCode;

    dest.ROIbox = source.ROIbox;
    dest.enforcelimits(dest.ROIbox);
    dest.activeROI = source.activeROI;
    if (dest.activeROI) {
      dest.Limits = source.Limits;  
      dest.enforcelimits(dest.Limits);
    } else {
      dest.setdefaultlimits();
    }
    dest.calc_no_voxels();
    
    dest.interpkernel = source.interpkernel;
    dest.p_interpmethod = source.p_interpmethod;
    dest.p_extrapmethod = source.p_extrapmethod;
    dest.padvalue = (D) source.padvalue;
    dest.splineorder = source.splineorder;
    dest.ep_valid = source.ep_valid;

    dest.displayMaximum=source.displayMaximum;
    dest.displayMinimum=source.displayMinimum;
    dest.setAuxFile(source.getAuxFile());
  }


  template <class S, class D>
  void copyconvert(const volume<S>& source, volume<D>& dest)
  {
    // set up basic size and data storage
    dest.reinitialize(source.xsize(),source.ysize(),source.zsize());
    // set up properties (except lazy ones)
    copybasicproperties(source,dest);
    // now copy across the data
    convertbuffer(source.Data,dest.Data,source.SizeBound);
    // doubly ensure that lazy properties are not valid
    dest.set_whole_cache_validity(false);
  }
  

  template <class S, class D>
  void copybasicproperties(const volume4D<S>& source, volume4D<D>& dest)
  {
    // set up properties (except lazy ones)
    dest.p_TR = source.p_TR;
    dest.ROIbox = source.ROIbox;
    dest.enforcelimits(dest.ROIbox);
    dest.activeROI = source.activeROI;
    if ((dest.activeROI) && (sameabssize(source,dest))) {
      dest.Limits = source.Limits;  
      dest.enforcelimits(dest.Limits);
    } else {
      dest.setdefaultlimits();
    }

    dest.p_interpmethod = source.p_interpmethod;
    dest.p_extrapmethod = source.p_extrapmethod;
    // Cannot convert the following between different types
    //dest.p_userextrap = source.p_userextrap;
    //dest.p_userinterp = source.p_userinterp;
    dest.p_padval = (D) source.p_padval;

    int toffset = dest.mint() - source.mint();
    for (int t=source.mint(); t<=source.maxt(); t++) {
      copybasicproperties(source[t],dest[Min(t + toffset,dest.maxt())]);
    }
  }

  template <class S, class D>
  void copybasicproperties(const volume4D<S>& source, volume<D>& dest)
  {
    copybasicproperties(source[0],dest);
  }

  template <class S, class D>
  void copybasicproperties(const volume<S>& source, volume4D<D>& dest)
  {
    // set up properties (except lazy ones)
    dest.setdefaultproperties();
    dest.setROIlimits(source.ROIlimits(0),source.ROIlimits(1),
		      source.ROIlimits(2),dest.ROIlimits(3),
		      source.ROIlimits(4),source.ROIlimits(5),
		      source.ROIlimits(6),dest.ROIlimits(7));
    if (source.usingROI()) dest.activateROI(); else dest.deactivateROI();
    if ((dest.usingROI()) && (dest.tsize()>=1) 
	&& (sameabssize(source,dest[0]))) {
      dest.setROIlimits(source.limits());  
    } else {
      dest.setdefaultlimits();
    }

    dest.setinterpolationmethod(source.getinterpolationmethod());
    dest.setextrapolationmethod(source.getextrapolationmethod());
    // Cannot convert the following between different types
    //dest.p_userextrap = source.p_userextrap;
    //dest.p_userinterp = source.p_userinterp;
    dest.setpadvalue((D) source.getpadvalue());

    for (int t=dest.mint(); t<=dest.maxt(); t++) {
      copybasicproperties(source,dest[t]);
    }
  }


  template <class S, class D>
  void copyconvert(const volume4D<S>& source, volume4D<D>& dest)
  {
    // set up basic size and data storage
    dest.reinitialize(source.xsize(),source.ysize(),source.zsize(),
		      source.tsize());

    copybasicproperties(source,dest);
    for (int t=0; t<source.tsize(); t++) {
      copyconvert(source[t],dest[t]);
    }

    // doubly ensure that lazy properties are not valid
    dest.set_whole_cache_validity(false);
  }
  


  template <class S>
  volume<S> operator+(S num, const volume<S>& vol)
    { return (vol + num); }

  template <class S>
  volume<S> operator-(S num, const volume<S>& vol)
  {
    volume<S> tmp = vol;
    tmp=num;
    tmp-=vol;
    return tmp;
  }

  template <class S>
  volume<S> operator*(S num, const volume<S>& vol)
    { return (vol * num); }

  template <class S>
  volume<S> operator/(S num, const volume<S>& vol)
  {
    volume<S> tmp = vol;
    tmp=num;
    tmp/=vol;
    return tmp;
  }

  template <class S>
  volume<S> operator-(const volume<S>& vol)
  {
    return(vol * (static_cast<S>(-1)));
  }

  template <class S>
  bool operator==(const volume<S>& v1,
		  const volume<S>& v2)
  {
    // Check relevant parts of header
    if (!samesize(v1,v2,true)) return(false);
    if (v1.sform_code() != v2.sform_code()) return(false);
    if (v1.sform_mat() != v2.sform_mat()) return(false);
    if (v1.qform_code() != v2.qform_code()) return(false);
    if (v1.qform_mat() != v2.qform_mat()) return(false);
    // Check data
    for (typename volume<S>::fast_const_iterator it1=v1.fbegin(), it_end=v1.fend(), it2=v2.fbegin(); it1 != it_end; ++it1, ++it2) {
      if ((*it1) != (*it2)) return(false);
    }

    return(true);
  }
  template <class S>
  bool operator!=(const volume<S>& v1, const volume<S>& v2) {return(!(v1==v2));}

  template <class S>
  volume4D<S> operator+(S num, const volume4D<S>& vol)
    { return (vol + num); }

  template <class S>
  volume4D<S> operator-(S num, const volume4D<S>& vol)
  {
    volume4D<S> tmp = vol;
    tmp=num;
    tmp-=vol;
    return tmp;
  }

  template <class S>
  volume4D<S> operator*(S num, const volume4D<S>& vol)
    { return (vol * num); }

  template <class S>
  volume4D<S> operator/(S num, const volume4D<S>& vol)
  {
    volume4D<S> tmp = vol;
    tmp=num;
    tmp/=vol;
    return tmp;
  }

  template <class S>
  volume4D<S> operator-(const volume4D<S>& vol)
  {
    return(vol * (static_cast<S>(-1)));
  }

  template <class S>
  volume4D<S> operator+(const volume<S>& v1, const volume4D<S>& v2)
    { return (v2 + v1); }

  template <class S>
  volume4D<S> operator-(const volume<S>& v1, const volume4D<S>& v2)
    { return (v2*(-1.0) + v1); }

  template <class S>
  volume4D<S> operator*(const volume<S>& v1, const volume4D<S>& v2)
    { return (v2 * v1); }

}  // end namespace

#endif








