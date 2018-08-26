//  
//  splinterpolator.h
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2008 University of Oxford 
//
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
//

#ifndef splinterpolator_h
#define splinterpolator_h

#include <vector>
#include <string>
#include <cmath>
#include "newmat.h"
#include "miscmaths/miscmaths.h"

namespace SPLINTERPOLATOR {

enum ExtrapolationType {Zeros, Constant, Mirror, Periodic};

class SplinterpolatorException: public std::exception
{
private:
  std::string m_msg;
public:
  SplinterpolatorException(const std::string& msg) throw(): m_msg(msg) {}

  virtual const char *what() const throw() {
    return string("Splinterpolator::" + m_msg).c_str();
  }

  ~SplinterpolatorException() throw() {}
};

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
// Class Splinterpolator:
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

template<class T>
class Splinterpolator
{
public:
  // Constructors
  Splinterpolator() : _valid(false), _own_coef(false), _coef(0), _cptr(0), _ndim(0) {}
  Splinterpolator(const T *data, const std::vector<unsigned int>& dim, const std::vector<ExtrapolationType>& et, unsigned int order=3, bool copy_low_order=true, double prec=1e-8) : _valid(false), _own_coef(false), _coef(0), _cptr(0), _ndim(0)
  {
    common_construction(data,dim,order,prec,et,copy_low_order);
  }
  Splinterpolator(const T *data, const std::vector<unsigned int>& dim, ExtrapolationType et=Zeros, unsigned int order=3, bool copy_low_order=true, double prec=1e-8) : _valid(false), _own_coef(false), _coef(0), _cptr(0), _ndim(0)
  {
    std::vector<ExtrapolationType>   ett(dim.size(),et);
    common_construction(data,dim,order,prec,ett,copy_low_order);
  }
  // Copy construction. May be removed in future
  Splinterpolator(const Splinterpolator<T>& src) : _valid(false), _own_coef(false), _coef(0), _cptr(0), _ndim(0) { assign(src); }

  // Destructor
  ~Splinterpolator() { if(_own_coef) delete [] _coef; }

  // Assignment. May be removed in future
  Splinterpolator& operator=(const Splinterpolator& src) { if(_own_coef) delete [] _coef; assign(src); return(*this); }

  // Set new data in Splinterpolator.
  void Set(const T *data, const std::vector<unsigned int>& dim, const std::vector<ExtrapolationType>& et, unsigned int order=3, bool copy_low_order=true, double prec=1e-8)
  {
    if (_own_coef) delete [] _coef;
    common_construction(data,dim,order,prec,et,copy_low_order);
  }
  void Set(const T *data, const std::vector<unsigned int>& dim, ExtrapolationType et, unsigned int order=3, bool copy_low_order=true, double prec=1e-8)
  {
    std::vector<ExtrapolationType>  vet(dim.size(),Zeros);
    Set(data,dim,vet,order,copy_low_order,prec);
  }

  // Return interpolated value
  T operator()(const std::vector<float>&  coord) const;
  T operator()(double x, double y=0, double z=0, double t=0) const
  {
   if (!_valid) throw SplinterpolatorException("operator(): Cannot interpolate un-initialized object");
   if (_ndim>4 || (t && _ndim<4) || (z && _ndim<3) || (y && _ndim<2)) throw SplinterpolatorException("operator(): input has wrong dimensionality");
    double coord[5] = {x,y,z,t,0.0};
    return(static_cast<T>(value_at(coord)));
  }

  // Return interpolated value along with first derivative in one direction (useful for distortion correction)
  T operator()(const std::vector<float>& coord, unsigned int dd, T *dval) const;
  T operator()(double x, double y, double z, unsigned int dd, T *dval) const; 
  T operator()(double x, double y, unsigned int dd, T *dval) const { return((*this)(x,y,0.0,dd,dval)); }
  T operator()(double x, T *dval) const { return((*this)(x,0.0,0.0,0,dval)); }

  // Return interpolated value along with selected derivatives
  T ValAndDerivs(const std::vector<float>& coord, const std::vector<unsigned int>& deriv, std::vector<T>& rderiv) const;
  T ValAndDerivs(const std::vector<float>& coord, std::vector<T>& rderiv) const
  {
    std::vector<unsigned int>   deriv(_ndim,1);
    return(ValAndDerivs(coord,deriv,rderiv));
  }
  T ValAndDerivs(double x, double y, double z, std::vector<T>& rderiv) const;

  // Return continous derivative at voxel centres (only works for order<1)
  T Deriv(const std::vector<unsigned int>& indx, unsigned int ddir) const;
  T Deriv1(const std::vector<unsigned int>& indx) const {return(Deriv(indx,0));}
  T Deriv2(const std::vector<unsigned int>& indx) const {return(Deriv(indx,1));}
  T Deriv3(const std::vector<unsigned int>& indx) const {return(Deriv(indx,2));}
  T Deriv4(const std::vector<unsigned int>& indx) const {return(Deriv(indx,3));}
  T Deriv5(const std::vector<unsigned int>& indx) const {return(Deriv(indx,4));}
  T DerivXYZ(unsigned int i, unsigned int j, unsigned int k, unsigned int dd) const;
  T DerivX(unsigned int i, unsigned int j, unsigned int k) const {return(DerivXYZ(i,j,k,0));}
  T DerivY(unsigned int i, unsigned int j, unsigned int k) const {return(DerivXYZ(i,j,k,1));}
  T DerivZ(unsigned int i, unsigned int j, unsigned int k) const {return(DerivXYZ(i,j,k,2));}
  void Grad3D(unsigned int i, unsigned int j, unsigned int k, T *xg, T *yg, T *zg) const;
  void Grad(const std::vector<unsigned int>& indx, std::vector<T>& grad) const;

  // Return continous addition (since previous voxel) of integral at voxel centres
  T IntX() const;
  T IntY() const;
  T IntZ() const;

  //
  // The "useful" functionality pretty much ends here.
  // Remaining functions are mainly for debugging/diagnostics.
  //
  unsigned int NDim() const { return(_ndim); }
  unsigned int Order() const { return(_order); }
  ExtrapolationType Extrapolation(unsigned int dim) const 
  { 
    if (dim >= _ndim) throw SplinterpolatorException("Extrapolation: Invalid dimension");
    return(_et[dim]);
  } 
  const std::vector<unsigned int>& Size() const { return(_dim); }
  unsigned int Size(unsigned int dim) const { if (dim > 4) return(0); else return(_dim[dim]);}
  T Coef(unsigned int x, unsigned int y=0, unsigned int z=0) const
  {
    std::vector<unsigned int> indx(3,0);
    indx[0] = x; indx[1] = y; indx[2] = z;
    return(Coef(indx));
  }
  T Coef(std::vector<unsigned int> indx) const;
  NEWMAT::ReturnMatrix CoefAsNewmatMatrix() const;
  NEWMAT::ReturnMatrix KernelAsNewmatMatrix(double sp=0.1, unsigned int deriv=0) const;
  
  //
  // Here we declare nested helper-class SplineColumn
  //
  class SplineColumn
  {
  public:
    // Constructor
    SplineColumn(unsigned int sz, unsigned int step) : _sz(sz), _step(step) { _col = new double[_sz]; }
    // Destructor
    ~SplineColumn() { delete [] _col; }
    // Extract a column from a volume
    void Get(const T  *dp)
    {
      for (unsigned int i=0; i<_sz; i++, dp+=_step) _col[i] = static_cast<double>(*dp);
    }
    // Insert column into volume
    void Set(T *dp) const
    {
      T   test = static_cast<T>(1.5);
      if (test == 1) {  // If T is not float or double
        for (unsigned int i=0; i<_sz; i++, dp+=_step) *dp = static_cast<T>(_col[i] + 0.5);   // Round to nearest integer
      }
      else {
        for (unsigned int i=0; i<_sz; i++, dp+=_step) *dp = static_cast<T>(_col[i]);
      }
    }
    // Deconvolve column
    void Deconv(unsigned int order, ExtrapolationType et, double prec);

  private:
    unsigned int  _sz;
    unsigned int  _step;
    double        *_col;

    unsigned int get_poles(unsigned int order, double *z, unsigned int *sf) const;
    double init_bwd_sweep(double z, double lv, ExtrapolationType et, double prec) const;
    double init_fwd_sweep(double z, ExtrapolationType et, double prec) const;

    SplineColumn(const SplineColumn& sc);              // Don't allow copy-construction
    SplineColumn& operator=(const SplineColumn& sc);   // Dont allow assignment
  };
  //
  // Here ends nested helper-class SplineColumn
  //

private:
  bool                            _valid;         // Decides if neccessary information has been set or not
  bool                            _own_coef;      // Decides if we "own" (have allocated) _coef
  T                               *_coef;         // Volume of spline coefficients
  const T                         *_cptr;         // Pointer to constant data. Used instead of _coef when we don't copy the data
  unsigned int                    _order;         // Order of splines
  unsigned int                    _ndim;          // # of non-singleton dimensions
  double                          _prec;          // Precision when dealing with edges
  std::vector<unsigned int>       _dim;           // Dimensions of data
  std::vector<ExtrapolationType>  _et;            // How to do extrapolation

  //
  // Private helper-functions
  //
  void common_construction(const T *data, const std::vector<unsigned int>& dim, unsigned int order, double prec, const std::vector<ExtrapolationType>& et, bool copy);
  void assign(const Splinterpolator<T>& src);
  bool calc_coef(const T *data, bool copy);
  void deconv_along(unsigned int dim);
  T coef(int *indx) const;
  const T* coef_ptr() const {if (_own_coef) return(_coef); else return(_cptr); }
  unsigned int indx2indx(int indx, unsigned int d) const;
  unsigned int indx2linear(int k, int l, int m) const;
  unsigned int add2linear(unsigned int lin, int j) const;
  double value_at(const double *coord) const;
  double value_and_derivatives_at(const double *coord, const unsigned int *deriv, double *dval) const; 
  void derivatives_at_i(const unsigned int *indx, const unsigned int *deriv, double *dval) const;
  unsigned int get_start_indicies(const double *coord, int *sinds) const;
  unsigned int get_start_indicies_at_i(const unsigned int *indx, int *sinds) const;
  unsigned int get_wgts(const double *coord, const int *sinds, double **wgts) const;
  unsigned int get_wgts_at_i(const unsigned int *indx, const int *sinds, double **wgts) const;
  unsigned int get_dwgts(const double *coord, const int *sinds, const unsigned int *deriv, double **dwgts) const;
  unsigned int get_dwgts_at_i(const unsigned int *indx, const int *sinds, const unsigned int *deriv, double **dwgts) const;
  double get_wgt(double x) const;
  double get_wgt_at_i(int i) const;
  double get_dwgt(double x) const;
  double get_dwgt_at_i(int i) const;
  void get_dwgt1(const double * const *wgts, const double * const *dwgts, const unsigned int *dd, unsigned int nd, 
                 unsigned int k, unsigned int l, unsigned int m, double wgt1, double *dwgt1) const;

  std::pair<double,double> range() const;
  bool should_be_zero(const double *coord) const;
  unsigned int n_nonzero(const unsigned int *vec) const;
  bool odd(unsigned int i) const {return(static_cast<bool>(i%2));}
  bool even(unsigned int i) const {return(!odd(i));}

  //
  // Disallowed member functions
  //
  // Splinterpolator(const Splinterpolator& s);             // Don't allow copy-construction
  // Splinterpolator& operator=(const Splinterpolator& s);  // Don't allow assignment
};

/////////////////////////////////////////////////////////////////////
//
// Here starts public member functions for Splinterpolator
//
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
//
// Returns interpolated value at location coord.
//
/////////////////////////////////////////////////////////////////////

template<class T>
T Splinterpolator<T>::operator()(const std::vector<float>& coord) const
{
  if (!_valid) throw SplinterpolatorException("operator(): Cannot interpolate un-initialized object");
  if (coord.size() != _ndim) throw SplinterpolatorException("operator(): coord has wrong length");
  double dcoord[5] = {0.0,0.0,0.0,0.0,0.0};
  for (unsigned int i=0; i<coord.size(); i++) dcoord[i] = coord[i];
  return(static_cast<T>(value_at(dcoord)));
}

/////////////////////////////////////////////////////////////////////
//
// Returns interpolated value and a single derivative at location coord.
// The derivative should be specified as the # of the dimension
// (starting at zero) that you want it along.
//
/////////////////////////////////////////////////////////////////////

template<class T>
T Splinterpolator<T>::operator()(const std::vector<float>& coord, unsigned int dd, T *dval) const
{
  if (!_valid) throw SplinterpolatorException("operator(): Cannot interpolate un-initialized object");
  if (coord.size() != _ndim) throw SplinterpolatorException("operator(): coord has wrong length");
  if (dd > (_ndim-1)) throw SplinterpolatorException("operator(): derivative specified for invalid direction");

  double          dcoord[5] = {0.0,0.0,0.0,0.0,0.0};
  for (unsigned int i=0; i<coord.size(); i++) dcoord[i] = coord[i];
  unsigned int    deriv[5] = {0,0,0,0,0};
  deriv[dd] = 1;
  double          ddval = 0.0;
  T               rval;
  rval = static_cast<T>(value_and_derivatives_at(dcoord,deriv,&ddval));
  *dval = static_cast<T>(ddval);
  
  return(rval);
}	   

/////////////////////////////////////////////////////////////////////
//
// Returns interpolated value and a single derivative at location
// given by x, y and . The derivative should be specified as the # 
// of the dimension (starting at zero) that you want it along.
//
/////////////////////////////////////////////////////////////////////

template<class T>
T Splinterpolator<T>::operator()(double x, double y, double z, unsigned int dd, T *dval) const
{
  if (!_valid) throw SplinterpolatorException("operator(): Cannot interpolate un-initialized object");
  if (_ndim>3 || (z && _ndim<3) || (y && _ndim<2)) throw SplinterpolatorException("operator(): input has wrong dimensionality");
  if (dd > (_ndim-1)) throw SplinterpolatorException("operator(): derivative specified for invalid direction");

  double          coord[5] = {x,y,z,0.0,0.0};
  unsigned int    deriv[5] = {0,0,0,0,0};
  deriv[dd] = 1;
  double          ddval = 0.0;
  T               rval;
  rval = static_cast<T>(value_and_derivatives_at(coord,deriv,&ddval));
  *dval = static_cast<T>(ddval);
  
  return(rval);
}	   

/////////////////////////////////////////////////////////////////////
//
// Returns interpolated value and selected (by deriv) derivatives 
// at location given by coord. The interpolated value is the return
// value and the derivatives are returned in rderiv. The input
// deriv should be an _ndim long vector where a 1 indicates that
// the derivative is required in that direction and a zero that it
// is not. 
//
/////////////////////////////////////////////////////////////////////

template<class T>
T Splinterpolator<T>::ValAndDerivs(const std::vector<float>& coord, const std::vector<unsigned int>& deriv, std::vector<T>& rderiv) const
{
  if (!_valid) throw SplinterpolatorException("ValAndDerivs: Cannot interpolate un-initialized object");
  if (coord.size() != _ndim || deriv.size() != _ndim) throw SplinterpolatorException("ValAndDerivs: input has wrong dimensionality");
  double        lcoord[5] = {0.0,0.0,0.0,0.0,0.0};
  unsigned int  lderiv[5] = {0,0,0,0,0};
  unsigned int  nd = 0;
  for (unsigned int i=0; i<coord.size(); i++) { lcoord[i] = coord[i]; nd += (lderiv[i]=(deriv[i])?1:0); }
  if (rderiv.size()!=nd) SplinterpolatorException("ValAndDerivs: mismatch between deriv and rderiv");
  double        dval[5];
  T rval = static_cast<T>(value_and_derivatives_at(lcoord,lderiv,dval));
  for (unsigned int i=0; i<nd; i++) rderiv[i] = static_cast<T>(dval[i]);
  return(rval);
}

/////////////////////////////////////////////////////////////////////
//
// Returns interpolated value and derivatives in the x, y and z
// directions at a location given by x, y and z. The interpolated 
// value is the return value and the derivatives are returned in rderiv. 
//
/////////////////////////////////////////////////////////////////////

template<class T>
T Splinterpolator<T>::ValAndDerivs(double x, double y, double z, std::vector<T>& rderiv) const
{
  if (!_valid) throw SplinterpolatorException("ValAndDerivs: Cannot interpolate un-initialized object");
  if (_ndim != 3 || rderiv.size() != _ndim) throw SplinterpolatorException("ValAndDerivs: input has wrong dimensionality");
  double        coord[5] = {x,y,z,0.0,0.0};
  unsigned int  deriv[5] = {1,1,1,0,0};
  double        dval[3];
  T rval = static_cast<T>(value_and_derivatives_at(coord,deriv,dval));
  for (unsigned int i=0; i<3; i++) rderiv[i] = static_cast<T>(dval[i]);
  return(rval);
}


/////////////////////////////////////////////////////////////////////
//
// Routine that returns a 3D gradient at an integer location.
//
/////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////
//
// Routine that returns a single derivative at an integer location.
//
/////////////////////////////////////////////////////////////////////

template<class T>
T Splinterpolator<T>::Deriv(const std::vector<unsigned int>& indx, unsigned int dd) const
{
  if (!_valid) throw SplinterpolatorException("Deriv: Cannot take derivative of un-initialized object");
  if (indx.size() != _ndim) SplinterpolatorException("Deriv: Input indx of wrong dimension");
  if (dd > (_ndim-1)) throw SplinterpolatorException("Deriv: derivative specified for invalid direction");
  double dval;
  unsigned int lindx[5] = {0,0,0,0,0};
  unsigned int deriv[5] = {0,0,0,0,0};
  for (unsigned int i=0; i<_ndim; i++) lindx[i]=indx[i];
  deriv[dd] = 1;
  derivatives_at_i(lindx,deriv,&dval);
  return(static_cast<T>(dval));
}

template<class T>
T Splinterpolator<T>::DerivXYZ(unsigned int i, unsigned int j, unsigned int k, unsigned int dd) const
{
  if (!_valid) throw SplinterpolatorException("DerivXYZ: Cannot take derivative of un-initialized object");
  if (_ndim!=3 || dd>2) throw SplinterpolatorException("DerivXYZ: Input has wrong dimensionality");
  double dval;
  unsigned int lindx[5] = {i,j,k,0,0};
  unsigned int deriv[5] = {0,0,0,0,0};
  deriv[dd] = 1;
  derivatives_at_i(lindx,deriv,&dval);
  return(static_cast<T>(dval));
}

template<class T>
void Splinterpolator<T>::Grad3D(unsigned int i, unsigned int j, unsigned int k, T *xg, T *yg, T *zg) const
{
  if (!_valid) throw SplinterpolatorException("Grad3D: Cannot take derivative of un-initialized object");
  if (_ndim != 3) SplinterpolatorException("Grad3D: Input of wrong dimension");
  unsigned int lindx[5] = {i,j,k,0,0};
  unsigned int deriv[5] = {1,1,1,0,0};
  double       dval[5] = {0.0,0.0,0.0,0.0,0.0};
  derivatives_at_i(lindx,deriv,dval);
  *xg=static_cast<T>(dval[0]); *yg=static_cast<T>(dval[1]); *zg=static_cast<T>(dval[2]);
  return; 
}

template<class T>
void Splinterpolator<T>::Grad(const std::vector<unsigned int>& indx, std::vector<T>& grad) const
{
  if (!_valid) throw SplinterpolatorException("Grad: Cannot take derivative of un-initialized object");
  if (indx.size() != _ndim || grad.size() != _ndim) SplinterpolatorException("Grad: Input indx or grad of wrong dimension");
  unsigned int lindx[5] = {0,0,0,0,0};
  unsigned int deriv[5] = {0,0,0,0,0};
  double       dval[5] = {0.0,0.0,0.0,0.0,0.0};
  for (unsigned int i=0; i<_ndim; i++) { lindx[i]=indx[i]; deriv[i]=1; }
  derivatives_at_i(lindx,deriv,dval);
  for (unsigned int i=0; i<_ndim; i++) grad[i] = static_cast<T>(dval[i]);
  return;
}

/////////////////////////////////////////////////////////////////////
//
// Returns the value of the coefficient given by indx (zero-offset)
//
/////////////////////////////////////////////////////////////////////

template<class T>
T Splinterpolator<T>::Coef(std::vector<unsigned int> indx) const
{
  if (!_valid) throw SplinterpolatorException("Coef: Cannot get coefficients for un-initialized object");
  if (!indx.size()) throw SplinterpolatorException("Coef: indx has zeros dimensions");
  if (indx.size() > 5) throw SplinterpolatorException("Coef: indx has more than 5 dimensions");
  for (unsigned int i=0; i<indx.size(); i++) if (indx[i] >= _dim[i]) throw SplinterpolatorException("Coef: indx out of range");

  unsigned int lindx=indx[indx.size()-1];
  for (int i=indx.size()-2; i>=0; i--) lindx = _dim[i]*lindx + indx[i];
  return(coef_ptr()[lindx]);
}

/////////////////////////////////////////////////////////////////////
//
// Returns the values of all coefficients as a Newmat matrix. If 
// _ndim==1 it will return a row-vector, if _ndim==2 it will return
// a matrix, if _ndim==3 it will return a tiled matrix where the n
// first rows (where n is the number of rows in one slice) pertain to 
// the first slice, the next n rows to the second slice, etc. And 
// correspondingly for 4- and 5-D.
//
/////////////////////////////////////////////////////////////////////

template<class T>
NEWMAT::ReturnMatrix Splinterpolator<T>::CoefAsNewmatMatrix() const
{
  if (!_valid) throw SplinterpolatorException("CoefAsNewmatMatrix: Cannot get coefficients for un-initialized object");
  NEWMAT::Matrix   mat(_dim[1]*_dim[2]*_dim[3]*_dim[4],_dim[0]);
  std::vector<unsigned int>  cindx(5,0);

  unsigned int r=0;
  for (cindx[4]=0; cindx[4]<_dim[4]; cindx[4]++) {
    for (cindx[3]=0; cindx[3]<_dim[3]; cindx[3]++) {
      for (cindx[2]=0; cindx[2]<_dim[2]; cindx[2]++) {
	for (cindx[1]=0; cindx[1]<_dim[1]; cindx[1]++, r++) {
	  for (cindx[0]=0; cindx[0]<_dim[0]; cindx[0]++) {
            mat.element(r,cindx[0]) = Coef(cindx);  
	  }
	}
      }
    }
  }
  mat.Release();
  return(mat);
}

/////////////////////////////////////////////////////////////////////
//
// Return the kernel matrix to verify its correctness.
//
/////////////////////////////////////////////////////////////////////

template<class T>
NEWMAT::ReturnMatrix Splinterpolator<T>::KernelAsNewmatMatrix(double sp,                // Distance (in ksp) between points 
                                                              unsigned int deriv) const // Derivative (only 0/1 implemented).
{
  if (!_valid) throw SplinterpolatorException("KernelAsNewmatMatrix: Cannot get kernel for un-initialized object");
  if (deriv > 1) throw SplinterpolatorException("KernelAsNewmatMatrix: only 1st derivatives implemented");

  std::pair<double,double>           rng = range();

  unsigned int i=0;
  for (double x=rng.first; x<=rng.second; x+=sp, i++) ; // Intentional
  NEWMAT::Matrix  kernel(i,2);
  for (double x=rng.first, i=0; x<=rng.second; x+=sp, i++) {
    kernel.element(i,0) = x;
    kernel.element(i,1) = (deriv) ? get_dwgt(x) : get_wgt(x);
  }
  kernel.Release();
  return(kernel);  
}
/////////////////////////////////////////////////////////////////////
//
// Here starts public member functions for SplineColumn
//
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
//
// This function implements the forward and backwards sweep
// as defined by equation 2.5 in Unsers paper:
//
// B-spline signal processing. II. Efficiency design and applications
//
/////////////////////////////////////////////////////////////////////

template<class T>
void Splinterpolator<T>::SplineColumn::Deconv(unsigned int order, ExtrapolationType et, double prec)
{
  double         z[3] = {0.0, 0.0, 0.0};     // Poles
  unsigned int   np = 0;                     // # of poles
  unsigned int   sf;                         // Scale-factor
  np = get_poles(order,z,&sf);

  for (unsigned int p=0; p<np; p++) {
    _col[0] = init_fwd_sweep(z[p],et,prec);
    double lv = _col[_sz-1];
    // Forward sweep
    double *ptr=&_col[1];
    for (unsigned int i=1; i<_sz; i++, ptr++) *ptr += z[p] * *(ptr-1);
    _col[_sz-1] = init_bwd_sweep(z[p],lv,et,prec);
    // Backward sweep
    ptr = &_col[_sz-2];
    for (int i=_sz-2; i>=0; i--, ptr--) *ptr = z[p]*(*(ptr+1) - *ptr);
  }
  double *ptr=_col;
  for (unsigned int i=0; i<_sz; i++, ptr++) *ptr *= sf; 
}

/////////////////////////////////////////////////////////////////////
//
// Here starts private member functions for Splinterpolator
//
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
//
// Returns the interpolated value at location given by coord. 
// coord must be a pointer to an array of indicies with _ndim
// values.
//
/////////////////////////////////////////////////////////////////////
/*
template<class T>
double Splinterpolator<T>::value_at(const double *coord) const
{
  if (should_be_zero(coord)) return(0.0);

  double       iwgt[8], jwgt[8], kwgt[8], lwgt[8], mwgt[8];
  double       *wgts[] = {iwgt, jwgt, kwgt, lwgt, mwgt};
  int          inds[5];
  unsigned int ni = 0;

  ni = get_start_indicies(coord,inds);
  get_wgts(coord,inds,wgts);

  double val=0.0;
  for (int m=0, me=(_ndim>4)?ni:1; m<me; m++) { 
    for (int l=0, le=(_ndim>3)?ni:1; l<le; l++) { 
      for (int k=0, ke=(_ndim>2)?ni:1; k<ke; k++) { 
        double wgt1 = wgts[4][m]*wgts[3][l]*wgts[2][k];
        for (int j=0, je=(_ndim>1)?ni:1; j<je; j++) { 
          double wgt2 = wgt1*wgts[1][j];
          for (int i=0; i<static_cast<int>(ni); i++) {
            int cindx[] = {inds[0]+i,inds[1]+j,inds[2]+k,inds[3]+l,inds[4]+m};
            val += coef(cindx)*wgts[0][i]*wgt2;
	  }
	}
      }
    }
  }
  return(val);
}
*/
template<class T>
double Splinterpolator<T>::value_at(const double *coord) const
{
  if (should_be_zero(coord)) return(0.0);

  double       iwgt[8], jwgt[8], kwgt[8], lwgt[8], mwgt[8];
  double       *wgts[] = {iwgt, jwgt, kwgt, lwgt, mwgt};
  int          inds[5];
  unsigned int ni = 0;
  const T      *cptr = coef_ptr();

  ni = get_start_indicies(coord,inds);
  get_wgts(coord,inds,wgts);

  double val=0.0;
  for (unsigned int m=0, me=(_ndim>4)?ni:1; m<me; m++) { 
    for (unsigned int l=0, le=(_ndim>3)?ni:1; l<le; l++) { 
      for (unsigned int k=0, ke=(_ndim>2)?ni:1; k<ke; k++) { 
        double wgt1 = wgts[4][m]*wgts[3][l]*wgts[2][k];
        unsigned int linear1 = indx2linear(inds[2]+k,inds[3]+l,inds[4]+m);
        for (unsigned int j=0, je=(_ndim>1)?ni:1; j<je; j++) { 
          double wgt2 = wgt1*wgts[1][j];
          int linear2 = add2linear(linear1,inds[1]+j);
          double *iiwgt=iwgt;
          for (unsigned int i=0; i<ni; i++, iiwgt++) {
	    val += cptr[linear2+indx2indx(inds[0]+i,0)]*(*iiwgt)*wgt2;
          }
	}
      }
    }
  }
  return(val);
}

/////////////////////////////////////////////////////////////////////
//
// Returns the interpolated value and selected derivatives at a 
// location given by coord. coord must be a pointer to an array 
// of voxel indicies with _ndim values. deriv must be a pointer
// to an _ndim long array of 0/1 specifying if the derivative is 
// requested in that direction or not.
//
/////////////////////////////////////////////////////////////////////

template<class T>
double Splinterpolator<T>::value_and_derivatives_at(const double       *coord,  
						    const unsigned int *deriv, 
						    double             *dval) 
const
{
  if (should_be_zero(coord)) { memset(dval,0,n_nonzero(deriv)*sizeof(double)); return(0.0); }

  double       iwgt[8], jwgt[8], kwgt[8], lwgt[8], mwgt[8];
  double       *wgts[] = {iwgt, jwgt, kwgt, lwgt, mwgt};
  double       diwgt[8], djwgt[8], dkwgt[8], dlwgt[8], dmwgt[8];
  double       *dwgts[] = {diwgt, djwgt, dkwgt, dlwgt, dmwgt};
  double       dwgt1[5];
  double       dwgt2[5];
  int          inds[5];
  unsigned int dd[5]; 
  unsigned int nd = 0;
  unsigned int ni = 0;
  const T      *cptr = coef_ptr();

  ni = get_start_indicies(coord,inds);
  get_wgts(coord,inds,wgts);
  get_dwgts(coord,inds,deriv,dwgts);
  for (unsigned int i=0; i<_ndim; i++) if (deriv[i]) { dd[nd] = i; dval[nd++] = 0.0; }

  double val=0.0;
  for (unsigned int m=0, me=(_ndim>4)?ni:1; m<me; m++) { 
    for (unsigned int l=0, le=(_ndim>3)?ni:1; l<le; l++) { 
      for (unsigned int k=0, ke=(_ndim>2)?ni:1; k<ke; k++) { 
        double wgt1 = wgts[4][m]*wgts[3][l]*wgts[2][k];
        get_dwgt1(wgts,dwgts,dd,nd,k,l,m,wgt1,dwgt1);
        unsigned int linear1 = indx2linear(inds[2]+k,inds[3]+l,inds[4]+m);
        for (unsigned int j=0, je=(_ndim>1)?ni:1; j<je; j++) { 
          double wgt2 = wgt1*wgts[1][j];
          for (unsigned int d=0; d<nd; d++) dwgt2[d] = (dd[d]==1) ? dwgt1[d]*dwgts[1][j] : dwgt1[d]*wgts[1][j];
          int linear2 = add2linear(linear1,inds[1]+j);
          double *iiwgt=iwgt;
          for (unsigned int i=0; i<ni; i++, iiwgt++) {
            double c = cptr[linear2+indx2indx(inds[0]+i,0)];
            val += c*(*iiwgt)*wgt2;
            for (unsigned int d=0; d<nd; d++) {
              double add = (dd[d]==0) ? c*diwgt[i]*dwgt2[d] : c*(*iiwgt)*dwgt2[d];
              dval[d] += add;
	    } 
	  }
	}
      }
    }
  }
  return(val);
}

template <class T>
void Splinterpolator<T>::derivatives_at_i(const unsigned int *indx,
                                          const unsigned int *deriv,
                                          double             *dval) 
const
{
  double       iwgt[8], jwgt[8], kwgt[8], lwgt[8], mwgt[8];
  double       *wgts[] = {iwgt, jwgt, kwgt, lwgt, mwgt};
  double       diwgt[8], djwgt[8], dkwgt[8], dlwgt[8], dmwgt[8];
  double       *dwgts[] = {diwgt, djwgt, dkwgt, dlwgt, dmwgt};
  double       dwgt1[5];
  double       dwgt2[5];
  int          inds[5];
  unsigned int dd[5]; 
  unsigned int nd = 0;
  unsigned int ni = 0;
  const T      *cptr = coef_ptr();

  ni = get_start_indicies_at_i(indx,inds);
  get_wgts_at_i(indx,inds,wgts);
  get_dwgts_at_i(indx,inds,deriv,dwgts);
  for (unsigned int i=0; i<_ndim; i++) if (deriv[i]) { dd[nd] = i; dval[nd++] = 0.0; }

  // double val=0.0;
  for (unsigned int m=0, me=(_ndim>4)?ni:1; m<me; m++) { 
    for (unsigned int l=0, le=(_ndim>3)?ni:1; l<le; l++) { 
      for (unsigned int k=0, ke=(_ndim>2)?ni:1; k<ke; k++) { 
        double wgt1 = wgts[4][m]*wgts[3][l]*wgts[2][k];
        get_dwgt1(wgts,dwgts,dd,nd,k,l,m,wgt1,dwgt1);
        unsigned int linear1 = indx2linear(inds[2]+k,inds[3]+l,inds[4]+m);
        for (unsigned int j=0, je=(_ndim>1)?ni:1; j<je; j++) { 
          // double wgt2 = wgt1*wgts[1][j];
          for (unsigned int d=0; d<nd; d++) dwgt2[d] = (dd[d]==1) ? dwgt1[d]*dwgts[1][j] : dwgt1[d]*wgts[1][j];
          int linear2 = add2linear(linear1,inds[1]+j);
          double *iiwgt=iwgt;
          for (unsigned int i=0; i<ni; i++, iiwgt++) {
            double c = cptr[linear2+indx2indx(inds[0]+i,0)];
            // val += c*(*iiwgt)*wgt2;
            for (unsigned int d=0; d<nd; d++) {
              double add = (dd[d]==0) ? c*diwgt[i]*dwgt2[d] : c*(*iiwgt)*dwgt2[d];
              dval[d] += add;
	    } 
	  }
	}
      }
    }
  }
  // return(val);
  return;
}
    
/////////////////////////////////////////////////////////////////////
//
// Returns (in sinds) the indicies of the first coefficient in all
// _ndim directions with a non-zero weight for the location given
// by coord. The caller is responsible for coord and sinds being
// valid pointers to arrays of 5 values.
// The return-value gives the total # of non-zero weights.
//
/////////////////////////////////////////////////////////////////////

template<class T>
unsigned int Splinterpolator<T>::get_start_indicies(const double *coord, int *sinds) const
{
  unsigned int ni = _order+1;
  
  if (odd(ni)) {
    for (unsigned int i=0; i<_ndim; i++) {
      sinds[i] = static_cast<int>(coord[i]+0.5) - ni/2;
    }
  }
  else {
    for (unsigned int i=0; i<_ndim; i++) {
      int ix = static_cast<int>(coord[i]+0.5);
      if (ix < coord[i]) sinds[i] = ix - (ni-1)/2;
      else sinds[i] = ix -ni/2;
    }
  }
  for (unsigned int i=_ndim; i<5; i++) sinds[i] = 0;

  return(ni);
}

// Does the same thing, but for integer (spot on voxel centre) index

template<class T>
unsigned int Splinterpolator<T>::get_start_indicies_at_i(const unsigned int *indx, int *sinds) const
{
  unsigned int ni = (odd(_order)) ? _order : _order+1;
  for (unsigned int i=0; i<_ndim; i++) {
    sinds[i] = indx[i] - (_order/2);
  }
  for (unsigned int i=_ndim; i<5; i++) sinds[i] = 0;
  return(ni);
}

/////////////////////////////////////////////////////////////////////
//
// Returns (in wgts) the weights for the coefficients given by sinds 
// for the location given by coord. 
//
/////////////////////////////////////////////////////////////////////

template<class T>
unsigned int Splinterpolator<T>::get_wgts(const double *coord, const int *sinds, double **wgts) const
{
  unsigned int ni = _order+1;

  for (unsigned int dim=0; dim<_ndim; dim++) {
    for (unsigned int i=0; i<ni; i++) {
      wgts[dim][i] = get_wgt(coord[dim]-(sinds[dim]+i));
    }
  }
  for (unsigned int dim=_ndim; dim<5; dim++) wgts[dim][0] = 1.0;

  return(ni);
}

// Same for integer (spot on voxel centre) index
template<class T>
unsigned int Splinterpolator<T>::get_wgts_at_i(const unsigned int *indx, const int *sinds, double **wgts) const
{
  unsigned int ni = (odd(_order)) ? _order : _order+1;
  
  for (unsigned int dim=0; dim<_ndim; dim++) {
    for (unsigned int i=0; i<ni; i++) {
      wgts[dim][i] = get_wgt_at_i(indx[dim]-(sinds[dim]+i));
    }
  }
  for (unsigned int dim=_ndim; dim<5; dim++) wgts[dim][0] = 1.0;

  return(ni);
}

template<class T>
unsigned int Splinterpolator<T>::get_dwgts(const double *coord, const int *sinds, const unsigned int *deriv, double **dwgts) const
{
  unsigned int ni = _order+1;

  for (unsigned int dim=0; dim<_ndim; dim++) {
    if (deriv[dim]) {
      switch (_order) {
      case 0:
	throw SplinterpolatorException("get_dwgts: invalid order spline");
        break;
      case 1:
        dwgts[dim][0] = -1; dwgts[dim][1] = 1;  // Not correct on original gridpoints
	break;
      case 2: case 3: case 4: case 5: case 6: case 7:
	for (unsigned int i=0; i<ni; i++) {
	  dwgts[dim][i] = get_dwgt(coord[dim]-(sinds[dim]+i));
        }
	break;
      default:
	throw SplinterpolatorException("get_dwgts: invalid order spline");
      }
    }
  }

  return(ni);
}

// Same for integer (spot on voxel centre) index
template<class T>
unsigned int Splinterpolator<T>::get_dwgts_at_i(const unsigned int *indx, const int *sinds, const unsigned int *deriv, double **dwgts) const
{
  unsigned int ni = (odd(_order)) ? _order : _order+1;

  for (unsigned int dim=0; dim<_ndim; dim++) {
    if (deriv[dim]) {
      switch (_order) {
      case 0: case 1:
	throw SplinterpolatorException("get_dwgts_at_i: invalid order spline");
        break;
      case 2: case 3: case 4: case 5: case 6: case 7:
	for (unsigned int i=0; i<ni; i++) {
	  dwgts[dim][i] = get_dwgt_at_i(indx[dim]-(sinds[dim]+i));
        }
	break;
      default:
	throw SplinterpolatorException("get_dwgts_at_i: invalid order spline");
      }
    }
  }

  return(ni);
}

/////////////////////////////////////////////////////////////////////
//
// Returns the weight for a spline at integer index i, where i is 
// relative to the centre index of the spline.
//
/////////////////////////////////////////////////////////////////////

template<class T>
double Splinterpolator<T>::get_wgt_at_i(int i) const
{
  double val = 0.0;
  int ai = std::abs(i);
  switch (_order) {
  case 0: case 1:
    val = (ai) ? 1.0 : 0.0;
    break;
  case 2:
    if (!ai) val = 0.75;
    else if (ai==1) val = 0.125;
    break;
  case 3:
    if (!ai) val = 0.666666666666667;
    else if (ai==1) val = 0.166666666666667;
    break;
  case 4:
    if (!ai) val = 0.598958333333333;
    else if (ai==1) val = 0.197916666666667;
    else if (ai==2) val = 0.002604166666667;
    break;
  case 5:
    if (!ai) val = 0.55;
    else if (ai==1) val = 0.216666666666667;
    else if (ai==2) val = 0.008333333333333;
    break;
  case 6:
    if (!ai) val = 0.511024305555556;
    else if (ai==1) val = 0.228797743055556;
    else if (ai==2) val = 0.015668402777779;
    else if (ai==3) val = 8.680555555555556e-05;
    break;
  case 7:
    if (!ai) val = 0.479365079365079;
    else if (ai==1) val = 0.236309523809524;
    else if (ai==2) val = 0.023809523809524;
    else if (ai==3) val = 1.984126984126984e-04;
    break;
  default:
    throw SplinterpolatorException("get_wgt_at_i: invalid order spline");
    break;
  }
  return(val);     
}

/////////////////////////////////////////////////////////////////////
//
// Returns the weight for the first derivative of a spline at integer 
// index i, where i is relative to the centre index of the spline.
//
/////////////////////////////////////////////////////////////////////

template<class T>
double Splinterpolator<T>::get_dwgt_at_i(int i) const
{
  double val = 0.0;
  int ai = std::abs(i);
  int sign = (ai) ? i/ai : 1;

  switch (_order) {
  case 0: case 1:
    throw SplinterpolatorException("get_dwgt: invalid order spline");
    break;
  case 2:
    if (!ai) val = 0.0;
    else if (ai==1) val = sign * (-0.5);
    break;
  case 3:
    if (!ai) val = 0.0;
    else if (ai==1) val = sign * (-0.5);
    break;
  case 4:
    if (!ai) val = 0.0;
    else if (ai==1) val = sign * (-0.458333333333333);
    else if (ai==2) val = sign * (-0.020833333333333);
    break;
  case 5:
    if (!ai) val = 0.0;
    else if (ai==1) val = sign * (-0.416666666666667);
    else if (ai==2) val = sign * (-0.041666666666667);
    break;
  case 6:
    if (!ai) val = 0.0;
    else if (ai==1) val = sign * (-0.376302083333333);
    else if (ai==2) val = sign * (-0.061458333333334);
    else if (ai==3) val = sign * (-2.604166666666667e-04);
    break;
  case 7:
    if (!ai) val = 0.0;
    else if (ai==1) val = sign * (-0.340277777777778);
    else if (ai==2) val = sign * (-0.077777777777778);
    else if (ai==3) val = sign * (-0.001388888888889);
    break;
  default:
    throw SplinterpolatorException("get_dwgt_at_i: invalid order spline");
    break;
  }
  return(val);     
}

/////////////////////////////////////////////////////////////////////
//
// Returns the weight for a spline at coordinate x, where x is relative
// to the centre of the spline.
//
/////////////////////////////////////////////////////////////////////

template<class T>
double Splinterpolator<T>::get_wgt(double x) const
{
  double val = 0.0;
  double ax = abs(x);  // Kernels all symmetric

  switch (_order) {
  case 0:
    if (ax < 0.5) val = 1.0;
    break;
  case 1:
    if (ax < 1) val = 1-ax;;
    break;
  case 2:
    if (ax < 0.5) val = 0.75-ax*ax;
    else if (ax < 1.5) val = 0.5*(1.5-ax)*(1.5-ax);
    break;
  case 3:
    if (ax < 1) val = 2.0/3.0 + 0.5*ax*ax*(ax-2);
    else if (ax < 2) { ax = 2-ax; val = (1.0/6.0)*(ax*ax*ax); }
    break;
  case 4:
    if (ax < 0.5) { ax *= ax; val = (115.0/192.0) + ax*((2.0*ax-5.0)/8.0); }
    else if (ax < 1.5) val = (55.0/96.0) + ax*(ax*(ax*((5.0-ax)/6.0) - 1.25) + 5.0/24.0);
    else if (ax < 2.5) { ax -= 2.5; ax *= ax; val = (1.0/24.0)*ax*ax; }
    break;
  case 5:
    if (ax < 1) { double xx = ax*ax; val = 0.55 + xx*(xx*((3.0-ax)/12.0) - 0.5); }
    else if (ax < 2) val = 0.425 + ax*(ax*(ax*(ax*((ax-9.0)/24.0) + 1.25) - 1.75) + 0.625);
    else if (ax < 3) { ax = 3-ax; double xx = ax*ax; val = (1.0/120.0)*ax*xx*xx; }
    break;
  case 6:
    if (ax < 0.5) { ax *= ax; val = (5887.0/11520.0) + ax*(ax*((21.0-4.0*ax)/144.0) -77.0/192.0); }
    else if (ax < 1.5) val = 7861.0/15360.0 + ax*(ax*(ax*(ax*(ax*((ax - 7.0)/48.0) + 0.328125) - 35.0/288.0) - 91.0/256.0) -7.0/768.0);
    else if (ax < 2.5) val = 1379.0/7680.0 + ax*(ax*(ax*(ax*(ax*((14.0-ax)/120.0) - 0.65625) + 133.0/72.0) - 2.5703125) + 1267.0/960.0);
    else if (ax < 3.5) { ax -= 3.5; ax *= ax*ax; val = (1.0/720.0) * ax*ax; }
    break;
  case 7:
    if (ax < 1) { double xx = ax*ax; val = 151.0/315.0 + xx*(xx*(xx*((ax-4.0)/144.0) + 1.0/9.0) - 1.0/3.0); }
    else if (ax < 2) val = 103.0/210.0 + ax*(ax*(ax*(ax*(ax*(ax*((12.0-ax)/240.0) -7.0/30.0) + 0.5) - 7.0/18.0) - 0.1) -7.0/90.0);
    else if (ax < 3) val = ax*(ax*(ax*(ax*(ax*(ax*((ax-20.0)/720.0) + 7.0/30.0) - 19.0/18.0) + 49.0/18.0) - 23.0/6.0) + 217.0/90.0) - 139.0/630.0;
    else if (ax < 4) { ax = 4-ax; double xxx=ax*ax*ax; val = (1.0/5040.0)*ax*xxx*xxx; }
    break;
  default:
    throw SplinterpolatorException("get_wgt: invalid order spline");
    break;
  }

  return(val);
}

/////////////////////////////////////////////////////////////////////
//
// Returns the weight for the first derivative of a spline at 
// coordinate x, where x is relative to the centre of the spline.
//
/////////////////////////////////////////////////////////////////////

template<class T>
double Splinterpolator<T>::get_dwgt(double x) const
{
  double val = 0.0;
  double ax = abs(x);                               // Kernels all anti-symmetric
  int    sign = (ax) ? static_cast<int>(x/ax) : 1;  // Arbitrary choice for when x=0

  switch (_order) {
  case 0: case 1:
    throw SplinterpolatorException("get_dwgt: invalid order spline");
    break;
  case 2:
    if (ax < 0.5) val = sign * -2.0*ax;
    else if (ax < 1.5) val = sign * (-1.5 + ax);
    break;
  case 3:
    if (ax < 1) val = sign * (1.5*ax*ax - 2.0*ax);
    else if (ax < 2) { ax = 2-ax; val = sign * -0.5*ax*ax; }
    break;
  case 4:
    if (ax < 0.5) val = sign * (ax*ax*ax - 1.25*ax);
    else if (ax < 1.5) val = sign * (5.0/24.0 - ax*(2.5 - ax*(2.5 - (2.0/3.0)*ax)));
    else if (ax < 2.5) { ax -= 2.5; val = sign * (1.0/6.0)*ax*ax*ax; }
    break;
  case 5:
    if (ax < 1) val = sign * ax*(ax*(ax*(1-(5.0/12.0)*ax)) - 1);
    else if (ax < 2) val = sign * (0.625 - ax*(3.5 - ax*(3.75 - ax*(1.5 - (5.0/24.0)*ax))));
    else if (ax < 3) { ax -= 3; ax = ax*ax; val = sign * (-1.0/24.0)*ax*ax; }
    break;
  case 6:
    if (ax < 0.5) { double xx = ax*ax; val = sign * ax*(xx*((7.0/12) - (1.0/6.0)*xx) - (77.0/96.0)); }
    else if (ax < 1.5) {double xx = ax*ax; val = sign * (ax*(xx*(0.1250*xx + 1.3125) - 0.7109375) - xx*((35.0/48.0)*xx + (35.0/96.0)) - (7.0/768.0)); }
    else if (ax < 2.5) { double xx = ax*ax; val = sign * ((1267.0/960.0) - ax*(xx*(0.05*xx + (21.0/8.0)) + (329.0/64.0)) + xx*((7.0/12.0)*xx + (133.0/24.0))); }
    else if (ax < 3.5) { ax -= 3.5; double xx = ax*ax; val = sign * (1.0/120.0)*xx*xx*ax; }
    break;
  case 7:
    if (ax < 1) { double xx = ax*ax; val = sign * ax*(xx*(xx*((7.0/144.0)*ax - (1.0/6.0)) + 4.0/9.0) - 2.0/3.0); }
    else if (ax < 2) { double xx = ax*ax; val = sign * (ax*(xx*(xx*0.3 + 2.0) - 0.2) - xx*(xx*(xx*(7.0/240.0) + (7.0/6.0)) + (7.0/6.0)) - (7.0/90.0)); }
    else if (ax < 3) { double xx = ax*ax; val = sign * (1.0/720.0)*(xx - 4.0*ax + 2.0)*(7.0*xx*xx - 92.0*xx*ax + 458.0*xx - 1024.0*ax + 868.0); }
    else if (ax < 4) { ax = 4-ax; ax = ax*ax*ax; val = sign * (-1.0/720.0)*ax*ax; }
    break;
  default:
    throw SplinterpolatorException("get_dwgt: invalid order spline");
    break;
  }

  return(val);
}

template<class T>
inline void Splinterpolator<T>::get_dwgt1(const double * const *wgts, const double * const *dwgts, 
                                          const unsigned int *dd, unsigned int nd, unsigned int k, 
                                          unsigned int l, unsigned int m, double wgt1, double *dwgt1) const
{
  for (unsigned int i=0; i<nd; i++) {
    switch (dd[i]) {
    case 2:
      dwgt1[i] = wgts[4][m] * wgts[3][l] * dwgts[2][k];
      break;
    case 3:
      dwgt1[i] = wgts[4][m] * dwgts[3][l] * wgts[2][k];
      break;
    case 4:
      dwgt1[i] = dwgts[4][m] * wgts[3][l] * wgts[2][k];
      break;
    default:
      dwgt1[i] = wgt1;
      break;
    }
  }
}

template<class T>
inline std::pair<double,double> Splinterpolator<T>::range() const
{
  std::pair<double,double>  rng(0.0,0.0);
  rng.second = static_cast<double>(_order+1.0)/2.0;
  rng.first = - rng.second;
  return(rng);
}
/////////////////////////////////////////////////////////////////////
//
// Returns the value of the coefficient indexed by indx. Unlike the
// public Coef() this routine allows indexing outside the valid
// volume, returning values that are dependent on the extrapolation
// model when these are encountered.
//
// N.B. May change value of input index N.B.
//
/////////////////////////////////////////////////////////////////////

template<class T>
inline unsigned int Splinterpolator<T>::indx2indx(int indx, unsigned int d) const
{
  if (d > (_ndim-1)) return(0);

  if (indx < 0) {
    switch (_et[d]) {
    case Constant:
      return(0);
      break;
    case Zeros: case Mirror:
      return(1-indx);
      break;
    case Periodic:
      return(_dim[d]+indx);
      break;
    default:
      break;
    }
  }
  else if (indx >= static_cast<int>(_dim[d])) {
    switch (_et[d]) {
    case Constant:
      return(_dim[d]-1);
      break;
    case Zeros: case Mirror:
      return(2*_dim[d]-indx-1);
      break;
    case Periodic:
      return(indx-_dim[d]);
      break;
    default:
      break;
    }
  }
  return(indx);
}

template<class T>
unsigned int Splinterpolator<T>::indx2linear(int k, int l, int m) const
{
  if (_ndim < 3) return(0);

  int lindx = 0;
  if (_ndim>4) lindx = indx2indx(m,4);
  if (_ndim>3) lindx = _dim[3]*lindx + indx2indx(l,3);
  lindx = _dim[0]*_dim[1]*(_dim[2]*lindx + indx2indx(k,2));

  return(lindx);
}

template<class T>
inline unsigned int Splinterpolator<T>::add2linear(unsigned int lin, int j) const
{
  if (_ndim < 2) return(lin);
  else return(lin + _dim[0]*indx2indx(j,1));
}

template<class T>
T Splinterpolator<T>::coef(int *indx) const
{
  // First fix any outside-volume indicies
  for (unsigned int i=0; i<_ndim; i++) {
    if (indx[i] < 0) {
      switch (_et[i]) {
      case Zeros:
	return(static_cast<T>(0));
	break;
      case Constant:
	indx[i] = 0;
	break;
      case Mirror:
	indx[i] = 1-indx[i];
	break;
      case Periodic:
	indx[i] = _dim[i]+indx[i];
	break;
      default:
	break;
      }
    }
    else if (indx[i] >= static_cast<int>(_dim[i])) {
      switch (_et[i]) {
      case Zeros:
	return(static_cast<T>(0));
	break;
      case Constant:
	indx[i] = _dim[i]-1;
	break;
      case Mirror:
	indx[i] = 2*_dim[i]-indx[i]-1;
	break;
      case Periodic:
	indx[i] = indx[i]-_dim[i];
	break;
      default:
	break;
      }
    }
  }
  // Now make linear index
  unsigned int lindx=indx[_ndim-1];
  for (int i=_ndim-2; i>=0; i--) lindx = _dim[i]*lindx + indx[i];
  return(coef_ptr()[lindx]);
}

template<class T>
bool Splinterpolator<T>::should_be_zero(const double *coord) const
{
  for (unsigned int i=0; i<_ndim; i++) {
    if (_et[i] == Zeros && (coord[i] < 0 || coord[i] > (_dim[i]-1))) return(true);
  }
  return(false);
}

template<class T>
unsigned int Splinterpolator<T>::n_nonzero(const unsigned int *vec) const
{
  unsigned int n=0;
  for (unsigned int i=0; i<_ndim; i++) if (vec[i]) n++;
  return(n);
}

/////////////////////////////////////////////////////////////////////
//
// Takes care of the "common" tasks when constructing a 
// Splinterpolator object. Called by constructors and by .Set()
//
/////////////////////////////////////////////////////////////////////

template<class T>
void Splinterpolator<T>::common_construction(const T *data, const std::vector<unsigned int>& dim, unsigned int order, double prec, const std::vector<ExtrapolationType>& et, bool copy)
{
  if (!dim.size()) throw SplinterpolatorException("common_construction: data has zeros dimensions");
  if (dim.size() > 5) throw SplinterpolatorException("common_construction: data cannot have more than 5 dimensions");
  if (dim.size() != et.size()) throw SplinterpolatorException("common_construction: dim and et must have the same size");
  for (unsigned int i=0; i<dim.size(); i++) if (!dim[i]) throw SplinterpolatorException("common_construction: data cannot have zeros size in any direction");
  if (order > 7) throw SplinterpolatorException("common_construction: spline order must be lesst than 7");
  if (!data) throw SplinterpolatorException("common_construction: zero data pointer");
  
  _order = order;
  _prec = prec;
  _et = et;
  _dim.resize(5);
  _ndim = dim.size();
  for (unsigned int i=0; i<5; i++) _dim[i]  = (i < dim.size()) ? dim[i] : 1;

  _own_coef = calc_coef(data,copy);

  _valid = true;
}

/////////////////////////////////////////////////////////////////////
//
// Takes care of the "common" tasks when copy-constructing
// and when assigning.
//
/////////////////////////////////////////////////////////////////////

template<class T>
void Splinterpolator<T>::assign(const Splinterpolator<T>& src)
{
  _valid = src._valid;
  _own_coef = src._own_coef;
  _cptr = src._cptr;
  _order = src._order;
  _ndim = src._ndim;
  _prec = src._prec;
  _dim = src._dim;
  _et = src._et;
  
  if (_own_coef) { // If we need to do a deep copy
    unsigned int ts = 1;
    for (unsigned int i=0; i<_ndim; i++) ts *= _dim[i];
    _coef = new T[ts];
    memcpy(_coef,src._coef,ts*sizeof(T));
  }
}

/////////////////////////////////////////////////////////////////////
//
// Performs deconvolution, converting signal to spline coefficients.
//
/////////////////////////////////////////////////////////////////////

template<class T>
bool Splinterpolator<T>::calc_coef(const T *data, bool copy)
{
  if (_order < 2 && !copy) { _cptr = data; return(false); }

  // Allocate memory and put the original data into _coef
  //
  unsigned int ts=1;
  for (unsigned int i=0; i<_dim.size(); i++) ts *= _dim[i];
  _coef = new T[ts];
  memcpy(_coef,data,ts*sizeof(T));

  if (_order < 2) return(true);  // If nearest neighbour or linear, that's all we need 

  // Loop over all non-singleton dimensions and deconvolve along them
  //
  std::vector<unsigned int>  tdim(_dim.size()-1,0); 
  for (unsigned int cdir=0; cdir<_dim.size(); cdir++) {
    if (_dim[cdir] > 1) deconv_along(cdir);
  }

  return(true);
}

/////////////////////////////////////////////////////////////////////
//
// Performs deconvolution along one of the dimensions, visiting
// all points along the other dimensions.
//
/////////////////////////////////////////////////////////////////////

template<class T>
void Splinterpolator<T>::deconv_along(unsigned int dim)
{
  // Set up to reflect "missing" dimension
  //
  std::vector<unsigned int>   rdim(4,1);     // Sizes along remaining dimensions
  std::vector<unsigned int>   rstep(4,1);    // Step-sizes (in "volume") of remaining dimensions
  unsigned int                mdim = 1;      // Size along "missing" dimension
  unsigned int                mstep = 1;     // Step-size along "missing" dimension
  for (unsigned int i=0, j=0, ss=1; i<5; i++) {
    if (i == dim) {  // If it is our "missing" dimension
      mdim = _dim[i];
      mstep = ss;
    }
    else {
      rdim[j] = _dim[i];
      rstep[j++] = ss;
    }
    ss *= _dim[i];
  }

  SplineColumn  col(mdim,mstep);          // Column helps us do the job
  
  for (unsigned int l=0; l<rdim[3]; l++) {
    for (unsigned int k=0; k<rdim[2]; k++) {
      for (unsigned int j=0; j<rdim[1]; j++) {
        T *dp = _coef + l*rstep[3] + k*rstep[2] + j*rstep[1];
	for (unsigned int i=0; i<rdim[0]; i++, dp+=rstep[0]) {
	  col.Get(dp);                          // Extract a column from the volume
          col.Deconv(_order,_et[dim],_prec);    // Deconvolve it
          col.Set(dp);                          // Put back the deconvolved column
	}
      }
    }
  }
  return;
}

/////////////////////////////////////////////////////////////////////
//
// Here starts private member functions for SplineColumn
//
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
//
// This function returns the poles and scale-factors for splines
// of order 2--7. The values correspond to those found in
// table 1 in Unsers 1993 paper: 
// B-spline signal processing. II. Efficiency design and applications
// 
// The actual values have been taken from
// http://bigwww.epfl.ch/thevenaz/interpolation/coeff.c
//
/////////////////////////////////////////////////////////////////////

template<class T>
unsigned int Splinterpolator<T>::SplineColumn::get_poles(unsigned int order, double *z, unsigned int *sf) const
{
  unsigned int   np = 0;                     // # of poles

  switch (order) {
  case 2:
    np = 1;
    z[0] = 2.0*sqrt(2.0) - 3.0;
    *sf = 8;
    break;
  case 3:
    np = 1;
    z[0] = sqrt(3.0) - 2.0;
    *sf = 6;
    break;
  case 4:
    np = 2;
    z[0] = sqrt(664.0 - sqrt(438976.0)) + sqrt(304.0) - 19.0;
    z[1] = sqrt(664.0 + sqrt(438976.0)) - sqrt(304.0) - 19.0;
    *sf = 384;
    break;
  case 5:
    np = 2;
    z[0] = sqrt(135.0 / 2.0 - sqrt(17745.0 / 4.0)) + sqrt(105.0 / 4.0) - 13.0 / 2.0;
    z[1] = sqrt(135.0 / 2.0 + sqrt(17745.0 / 4.0)) - sqrt(105.0 / 4.0) - 13.0 / 2.0;
    *sf = 120;
    break;
  case 6:
    np = 3;
    z[0] = -0.48829458930304475513011803888378906211227916123938;
    z[1] = -0.081679271076237512597937765737059080653379610398148;
    z[2] = -0.0014141518083258177510872439765585925278641690553467;
    *sf = 46080;
    break;
  case 7:
    np = 3;
    z[0] = -0.53528043079643816554240378168164607183392315234269;
    z[1] = -0.12255461519232669051527226435935734360548654942730;
    z[2] = -0.0091486948096082769285930216516478534156925639545994;
    *sf = 5040;
    break;
  default:
    throw SplinterpolatorException("SplineColumn::get_poles: invalid order of spline");
    break;
  }
  return(np);
}

/////////////////////////////////////////////////////////////////////
//
// Initialises the first value for the forward sweep. The initialisation
// will always be an approximation (this is where the "infinite" in IIR
// breaks down) so the value will be calculated to a predefined precision.
//
/////////////////////////////////////////////////////////////////////

template<class T>
double Splinterpolator<T>::SplineColumn::init_fwd_sweep(double z, ExtrapolationType et, double prec) const
{
  //
  // Move logs away from here after debugging
  //
  unsigned int n = static_cast<unsigned int>((log(prec)/log(abs(z))) + 1.5);
  n = (n > _sz) ? _sz : n;

  double iv = _col[0];
  if (et == Periodic) {
    double *ptr=&_col[_sz-1];
    double z2i=z;
    for (unsigned int i=1; i<n; i++, ptr--, z2i*=z) iv += z2i * *ptr;
  }
  else {
    double *ptr=&_col[1];
    double z2i=z;
    for (unsigned int i=1; i<n; i++, ptr++, z2i*=z) iv += z2i * *ptr;
  }
  return(iv); 
}

/////////////////////////////////////////////////////////////////////
//
// Initialises the first value for the backward sweep. The initialisation
// will always be an approximation (this is where the "infinite" in IIR
// breaks down) so the value will be calculated to a predefined precision.
//
/////////////////////////////////////////////////////////////////////

template<class T>
double Splinterpolator<T>::SplineColumn::init_bwd_sweep(double z, double lv, ExtrapolationType et, double prec) const
{
  double iv = 0.0;
  if (et == Periodic) {
    unsigned int n = static_cast<unsigned int>((log(prec)/log(abs(z))) + 1.5);
    n = (n > _sz) ? _sz : n;
    iv = z * _col[_sz-1];
    double z2i = z*z;
    double *ptr=_col;
    for (unsigned int i=1; i<n; i++, ptr++, z2i*=z) {
      iv += z2i * *ptr;
    }
    iv /= (z2i-1.0);
  }
  else {
    iv = -z/(1.0-z*z) * (2.0*_col[_sz-1] - lv);
  }
  return(iv);
}

} // End namespace SPLINTERPOLATOR

#endif // End #ifndef splinterpolator.h
