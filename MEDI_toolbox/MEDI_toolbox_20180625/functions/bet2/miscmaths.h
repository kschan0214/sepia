/*  miscmaths.h

    Mark Jenkinson & Mark Woolrich & Christian Beckmann & Tim Behrens, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

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

// Miscellaneous maths functions




#if !defined(__miscmaths_h)
#define __miscmaths_h

#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <assert.h>
#include <fstream>
#include <sstream>
#include <string>
#include <iterator>
#include "fslio/fslio.h"
//#include "config.h"
#include "newmatap.h"
#include "kernel.h"

//#pragma interface

using namespace NEWMAT;
using namespace std;

namespace MISCMATHS {

#ifndef M_PI
#define M_PI            3.14159265358979323846
#endif

#define OUT(t) cout<<#t "="<<t<<endl;
#define LOGOUT(t) LogSingleton::getInstance().str()<<#t "="<<t<<endl;

  // IO/string stuff
  template <class T> string num2str(T n, int width=-1);

#if defined(_MSC_VER) && (_MSC_VER < 1300)
  template <class T> string num2str(T n) { return num2str(n -1); }
#endif

  string size(const Matrix& mat);
  bool isNumber(const string& str);
  ReturnMatrix read_ascii_matrix(const string& filename, int nrows, int ncols);
  ReturnMatrix read_ascii_matrix(int nrows, int ncols, const string& filename);
  ReturnMatrix read_ascii_matrix(const string& filename);
  ReturnMatrix read_vest(string p_fname);
  int read_binary_matrix(Matrix& mres, const string& filename);
  ReturnMatrix read_binary_matrix(const string& filename);
  ReturnMatrix read_matrix(const string& filename);

  int write_ascii_matrix(const Matrix& mat, const string& filename, 
			 int precision=-1);
  int write_ascii_matrix(const string& filename, const Matrix& mat, 
			 int precision=-1);
  int write_vest(const Matrix& x, string p_fname, int precision=-1);
  int write_vest(string p_fname, const Matrix& x, int precision=-1);
  int write_binary_matrix(const Matrix& mat, const string& filename);

  // more basic IO calls
  string skip_alpha(ifstream& fs);
  ReturnMatrix read_ascii_matrix(ifstream& fs, int nrows, int ncols);
  ReturnMatrix read_ascii_matrix(int nrows, int ncols, ifstream& fs);
  ReturnMatrix read_ascii_matrix(ifstream& fs);
  int read_binary_matrix(Matrix& mres, ifstream& fs);
  ReturnMatrix read_binary_matrix(ifstream& fs);
  int write_ascii_matrix(const Matrix& mat, ofstream& fs, int precision=-1);
  int write_ascii_matrix(ofstream& fs, const Matrix& mat, int precision=-1);
  int write_binary_matrix(const Matrix& mat, ofstream& fs);

  // General maths

  int round(int x); 
  int round(float x);
  int round(double x);
  double rounddouble(double x);

  inline int sign(int x){ if (x>0) return 1; else { if (x<0) return -1; else return 0; } }
  inline int sign(float x){ if (x>0) return 1; else { if (x<0) return -1; else return 0; } }
  inline int sign(double x){ if (x>0) return 1; else { if (x<0) return -1; else return 0; } }
  
  inline double pow(double x, float y) { return std::pow(x,(double) y); }
  inline double pow(float x, double y) { return std::pow((double) x,y); }
  inline double pow(double x, int y) { return std::pow(x,(double) y); }
  inline float pow(float x, int y) { return std::pow(x,(float) y); }
  inline double pow(int x, double y) { return std::pow((double)x, y); }
  inline float pow(int x, float y) { return std::pow((float)x, y); }

  inline double sqrt(int x) { return std::sqrt((double) x); }
  inline double log(int x) { return std::log((double) x); }

  float Sinc(const float x);
  double Sinc(const double x);

  int periodicclamp(int x, int x1, int x2);

  template<class S, class T> 
   inline T Min(const S &a, const T &b) { if (a<b) return (T) a; else return b; }

  template<class S, class T>
   inline T Max(const S &a, const T &b) { if (a>b) return (T) a; else return b; }

  template<class T>
   inline T Sqr(const T& x) { return x*x; } 

  ColumnVector cross(const ColumnVector& a, const ColumnVector& b);
  ColumnVector cross(const Real *a, const Real *b);

  inline float dot(const ColumnVector& a, const ColumnVector& b)
    { return Sum(SP(a,b)); }

  double norm2(const ColumnVector& x);
  double norm2sq(double a, double b, double c);
  float norm2sq(float a, float b, float c);

  ColumnVector seq(const int num);

  int diag(Matrix& m, const float diagvals[]);
  int diag(Matrix& m, const ColumnVector& diagvals);
  int diag(DiagonalMatrix& m, const ColumnVector& diagvals);
  ReturnMatrix diag(const Matrix& mat);
  ReturnMatrix pinv(const Matrix& mat);
  int rank(const Matrix& X);
  ReturnMatrix sqrtaff(const Matrix& mat);

  void reshape(Matrix& r, const Matrix& m, int nrows, int ncols);
  ReturnMatrix reshape(const Matrix& m, int nrows, int ncols);
  int addrow(Matrix& m, int ncols);

  int construct_rotmat_euler(const ColumnVector& params, int n, Matrix& aff);
  int construct_rotmat_euler(const ColumnVector& params, int n, Matrix& aff,
		   const ColumnVector& centre);
  int construct_rotmat_quat(const ColumnVector& params, int n, Matrix& aff);
  int construct_rotmat_quat(const ColumnVector& params, int n, Matrix& aff,
		   const ColumnVector& centre);
  int make_rot(const ColumnVector& angl, const ColumnVector& centre, 
	       Matrix& rot);

  int getrotaxis(ColumnVector& axis, const Matrix& rotmat);
  int rotmat2euler(ColumnVector& angles, const Matrix& rotmat);
  int rotmat2quat(ColumnVector& quaternion, const Matrix& rotmat);
  int decompose_aff(ColumnVector& params, const Matrix& affmat,
		    int (*rotmat2params)(ColumnVector& , const Matrix& ));
  int decompose_aff(ColumnVector& params, const Matrix& affmat,
		    const ColumnVector& centre,
		    int (*rotmat2params)(ColumnVector& , const Matrix& ));
  int compose_aff(const ColumnVector& params, int n, const ColumnVector& centre,
		  Matrix& aff, 
		  int (*params2rotmat)(const ColumnVector& , int , Matrix& ,
				       const ColumnVector& ) );
  float rms_deviation(const Matrix& affmat1, const Matrix& affmat2, 
		      const ColumnVector& centre, const float rmax); 
  float rms_deviation(const Matrix& affmat1, const Matrix& affmat2, 
		      const float rmax=80.0); 

  Matrix Mat44ToNewmat(mat44 m);
  mat44 NewmatToMat44(const Matrix& m);
  mat44 newmat_to_mat44(const Matrix& inmat);
  Matrix mat44_to_newmat(mat44 inmat);

  void get_axis_orientations(const Matrix& sform_mat, int sform_code,
			     const Matrix& qform_mat, int qform_code,
			     int& icode, int& jcode, int& kcode);
    
  // 1D lookup table with linear interpolation
  float interp1(const ColumnVector& x, const ColumnVector& y, float xi);

  float quantile(const ColumnVector& in, int which);
  float percentile(const ColumnVector& in, float p);
  inline float median(const ColumnVector& in){ return quantile(in,2);}
  inline float iqr(const ColumnVector &in) { return quantile(in,3) - quantile(in,1); }

  ReturnMatrix quantile(const Matrix& in, int which);
  ReturnMatrix percentile(const Matrix& in, float p);
  inline ReturnMatrix median(const Matrix& in){ return quantile(in,2);}
  inline ReturnMatrix iqr(const Matrix& in){ Matrix res = quantile(in,3) - quantile(in,1); res.Release(); return res;}

  void cart2sph(const ColumnVector& dir, float& th, float& ph);// cartesian to sperical polar coordinates
  void cart2sph(const Matrix& dir,ColumnVector& th,ColumnVector& ph);//ditto
  void cart2sph(const vector<ColumnVector>& dir,ColumnVector& th,ColumnVector& ph);// same but in a vector

  // geometry function
  inline float point_plane_distance(const ColumnVector& X,const ColumnVector& P){//plane defined by a,b,c,d with a^2+b^2+c^2=1
    return( dot(X,P.SubMatrix(1,3,1,1))+P(4) );
  }

  // returns the first P such that 2^P >= abs(N). 
  int nextpow2(int n);

  // Auto-correlation function estimate of columns of p_ts
  // gives unbiased estimate - scales the raw correlation by 1/(N-abs(lags))
  void xcorr(const Matrix& p_ts, Matrix& ret, int lag = 0, int p_zeropad = 0);
  ReturnMatrix xcorr(const Matrix& p_ts, int lag = 0, int p_zeropad = 0);

  // removes trend from columns of p_ts
  // if p_level==0 it just removes the mean
  // if p_level==1 it removes linear trend
  // if p_level==2 it removes quadratic trend
  void detrend(Matrix& p_ts, int p_level=1);

  ReturnMatrix zeros(const int dim1, const int dim2 = -1);
  ReturnMatrix ones(const int dim1, const int dim2 = -1);
  ReturnMatrix repmat(const Matrix& mat, const int rows = 1, const int cols = 1);
  ReturnMatrix dist2(const Matrix& mat1, const Matrix& mat2);
  ReturnMatrix abs(const Matrix& mat);
  ReturnMatrix sqrt(const Matrix& mat);
  ReturnMatrix sqrtm(const Matrix& mat);
  ReturnMatrix log(const Matrix& mat);
  ReturnMatrix exp(const Matrix& mat);
  ReturnMatrix expm(const Matrix& mat);
  ReturnMatrix tanh(const Matrix& mat);
  ReturnMatrix pow(const Matrix& mat, const double exp);
  ReturnMatrix sum(const Matrix& mat, const int dim = 1);
  ReturnMatrix mean(const Matrix& mat, const int dim = 1);
  ReturnMatrix var(const Matrix& mat, const int dim = 1);
  ReturnMatrix max(const Matrix& mat);
  ReturnMatrix max(const Matrix& mat,ColumnVector& index);
  ReturnMatrix min(const Matrix& mat);
  ReturnMatrix gt(const Matrix& mat1,const Matrix& mat2); 
  ReturnMatrix lt(const Matrix& mat1,const Matrix& mat2); 
  ReturnMatrix geqt(const Matrix& mat1,const Matrix& mat2);  
  ReturnMatrix geqt(const Matrix& mat1,const float a); 
  ReturnMatrix leqt(const Matrix& mat1,const Matrix& mat2); 
  ReturnMatrix eq(const Matrix& mat1,const Matrix& mat2); 
  ReturnMatrix neq(const Matrix& mat1,const Matrix& mat2); 
  ReturnMatrix SD(const Matrix& mat1,const Matrix& mat2); // Schur (element-wise) divide
  ReturnMatrix vox_to_vox(const ColumnVector& xyz1,const ColumnVector& dims1,const ColumnVector& dims2,const Matrix& xfm);
  ReturnMatrix mni_to_imgvox(const ColumnVector& mni,const ColumnVector& mni_origin,const Matrix& mni2img, const ColumnVector& img_dims);
  void remmean(const Matrix& mat, Matrix& demeanedmat, Matrix& Mean,  const int dim = 1);
  ReturnMatrix remmean(const Matrix& mat, const int dim = 1);
  ReturnMatrix stdev(const Matrix& mat, const int dim = 1);
  ReturnMatrix cov(const Matrix& mat, const int norm = 0);
  ReturnMatrix corrcoef(const Matrix& mat, const int norm = 0);
  void symm_orth(Matrix &Mat);
  void powerspectrum(const Matrix &Mat1, Matrix &Result, bool useLog);
  void element_mod_n(Matrix& Mat,double n); //represent each element in modulo n (useful for wrapping phases (n=2*M_PI))

  // matlab-like flip function
  ReturnMatrix flipud(const Matrix& mat);
  ReturnMatrix fliplr(const Matrix& mat);
  
  // ols
  // data is t x v
  // des is t x ev (design matrix)
  // tc is cons x ev (contrast matrix)
  // cope and varcope will be cons x v
  // but will be resized if they are wrong
  void ols(const Matrix& data,const Matrix& des,const Matrix& tc, Matrix& cope,Matrix& varcope);
  float ols_dof(const Matrix& des);
  

  // Conjugate Gradient methods to solve for x in:   A * x = b 
  // A must be symmetric and positive definite
  int conjgrad(ColumnVector& x, const Matrix& A, const ColumnVector& b, 
	       int maxit=3);
  // allow specification of reltol = relative tolerance of residual error
  //  (stops when error < reltol * initial error)
  int conjgrad(ColumnVector& x, const Matrix& A, const ColumnVector& b, 
	       int maxit, float reltol);

  float csevl(const float x, const ColumnVector& cs, const int n);
  float digamma(const float x);
  void glm_vb(const Matrix& X, const ColumnVector& Y, ColumnVector& B, SymmetricMatrix& ilambda_B, int niters=20);

  vector<float> ColumnVector2vector(const ColumnVector& col);
  
  ///////////////////////////////////////////////////////////////////////////
  // Uninteresting byte swapping functions
  void Swap_2bytes ( int n , void *ar ) ;
  void Swap_4bytes ( int n , void *ar ) ;
  void Swap_8bytes ( int n , void *ar ) ;
  void Swap_16bytes( int n , void *ar ) ;
  void Swap_Nbytes ( int n , int siz , void *ar ) ;

  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////

  // TEMPLATE DEFINITIONS //

  template<class t> ReturnMatrix vector2ColumnVector(const vector<t>& vec)
  {
    ColumnVector col(vec.size());
    
    for(unsigned int c = 0; c < vec.size(); c++)
      col(c+1) = vec[c];

    col.Release();
    return col;
  }
 
  template<class t> void write_vector(const string& fname, const vector<t>& vec)
  { 
    ofstream out;
    out.open(fname.c_str(), ios::out);
    copy(vec.begin(), vec.end(), ostream_iterator<t>(out, " "));
  }
  
  template<class t> void write_vector(const vector<t>& vec, const string& fname)
  {
    write_vector(fname,vec);
  }  
  
  template <class T>
  string num2str(T n, int width)
  {
    ostringstream os;
    if (width>0) {
      os.fill('0');
      os.width(width);
      os.setf(ios::internal, ios::adjustfield);
    }
    os << n;
    return os.str();
  }

}

#endif
