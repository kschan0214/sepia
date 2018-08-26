/*  miscmaths.cc

    Mark Jenkinson, Mark Woolrich, Christian Beckmann, Tim Behrens and Matthew Webster, FMRIB Image Analysis Group

    Copyright (C) 1999-2009 University of Oxford  */

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
#define NOMINMAX
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include "miscmaths.h"
#include "miscprob.h"
#include "stdlib.h"
#include "newmatio.h"

using namespace std;

namespace MISCMATHS {

  // The following lines are ignored by the current SGI compiler
  //  (version egcs-2.91.57)
  // A temporary fix of including the std:: in front of all abs() etc
  //  has been done for now
  using std::abs;
  using std::sqrt;
  using std::exp;
  using std::log;
  //  using std::pow;
  using std::atan2;


  string size(const Matrix& mat)
  {
    string str = num2str(mat.Nrows())+"*"+num2str(mat.Ncols());
        
    return str;
  }



  float Sinc(const float x) {
    if (fabs(x)<1e-9) {
      return 1-x*x*M_PI*M_PI/6.0;
    } else {
      return sin(M_PI*x)/(M_PI*x);
    }
  }

  double Sinc(const double x) {
    if (fabs(x)<1e-9) {
      return 1-x*x*M_PI*M_PI/6.0;
    } else {
      return sin(M_PI*x)/(M_PI*x);
    }
  }

  // General string/IO functions
  bool isNumber( const string& input)
  {
    if (input.size()==0) return false;
    char *pend;
    strtod(input.c_str(),&pend);
    if (*pend!='\0') return false;
    return true; 
  } 

  string skip_alpha(ifstream& fs) 
  {
    string cline;
    while (!fs.eof()) {
      streampos curpos = fs.tellg();
      getline(fs,cline);
      cline += " "; // force extra entry in parsing
      istringstream ss(cline.c_str());
      string firstToken="";
      ss >> firstToken; //Put first non-whitespace sequence into cc
      if (isNumber(firstToken)) {
	if (!fs.eof()) { fs.seekg(curpos); } else { fs.clear(); fs.seekg(0,ios::beg); }
	return cline;
      }
    }
    return "";
  }

  ReturnMatrix read_ascii_matrix(int nrows, int ncols, const string& filename)
  {
    return read_ascii_matrix(filename,nrows,ncols);
  }


  ReturnMatrix read_ascii_matrix(const string& filename, int nrows, int ncols)
  {
    Matrix mat(nrows,ncols);
    mat = 0.0;
  
    if ( filename.size()<1 ) return mat;
    ifstream fs(filename.c_str());
    if (!fs) { 
      cerr << "Could not open matrix file " << filename << endl;
      return mat;
    }
    mat = read_ascii_matrix(fs,nrows,ncols);
    fs.close();
    mat.Release();
    return mat;
  }

  ReturnMatrix read_ascii_matrix(int nrows, int ncols, ifstream& fs)
  {
    return read_ascii_matrix(fs, nrows, ncols);
  }

  ReturnMatrix read_ascii_matrix(ifstream& fs, int nrows, int ncols)
  {
    Matrix mat(nrows,ncols);
    mat = 0.0;
    string ss="";
  
    ss = skip_alpha(fs);
    for (int r=1; r<=nrows; r++) {
      for (int c=1; c<=ncols; c++) {
	if (!fs.eof()) {
	  fs >> ss;
	  while ( !isNumber(ss) && !fs.eof() ) {
	    fs >> ss;
	  }
	  mat(r,c) = atof(ss.c_str());
	}
      }
    }
    mat.Release();
    return mat;
  }

  
  ReturnMatrix read_ascii_matrix(const string& filename)
  {
    Matrix mat;
    if ( filename.size()<1 ) return mat;
    ifstream fs(filename.c_str());
    if (!fs) { 
      cerr << "Could not open matrix file " << filename << endl;
      mat.Release();
      return mat;
    }
    mat = read_ascii_matrix(fs);
    fs.close();
    mat.Release();
    return mat;
  }


  ReturnMatrix read_ascii_matrix(ifstream& fs)
  {
    int nRows(0), nColumns(0);
    string currentLine;
    // skip initial non-numeric lines
    //  and count the number of columns in the first numeric line
    currentLine = skip_alpha(fs);
    currentLine += " ";
    {
      istringstream ss(currentLine.c_str());
      string dummyToken="";
      while (!ss.eof()) {
	nColumns++;
	ss >> dummyToken;
      }
    }
    nColumns--;

    do {
      getline(fs,currentLine);
      currentLine += " "; // force extra entry in parsing
      istringstream ss(currentLine.c_str()); 
      string firstToken("");
      ss >> firstToken; //Put first non-whitespace sequence into cc
      if (!isNumber(firstToken)) break;  // stop processing when non-numeric line found
      nRows++;  // add new row to matrix
    } while (!fs.eof());
    
    // now know the size of matrix
    fs.clear();
    fs.seekg(0,ios::beg);    
    return read_ascii_matrix(fs,nRows,nColumns);

  }

#define BINFLAG 42

  ReturnMatrix read_binary_matrix(const string& filename)
  {
    Matrix mres;
    read_binary_matrix(mres,filename);
    mres.Release();
    return mres;
  }


  int read_binary_matrix(Matrix& mres, const string& filename)
  {
    if ( filename.size()<1 ) return 1;
    ifstream fs(filename.c_str(), ios::in | ios::binary);
    if (!fs) { 
      cerr << "Could not open matrix file " << filename << endl;
      return 2;
    }
    read_binary_matrix(mres,fs);
    fs.close();
    return 0;
  }

  ReturnMatrix read_binary_matrix(ifstream& fs)
  {
    Matrix mres;
    read_binary_matrix(mres,fs);
    mres.Release();
    return mres;
  }

  int read_binary_matrix(Matrix& mres, ifstream& fs)
  {
    bool swapbytes = false;
    unsigned int testval;
    // test for byte swapping
    fs.read((char*)&testval,sizeof(testval));
    if (testval!=BINFLAG) {
      swapbytes = true;
      Swap_Nbytes(1,sizeof(testval),&testval);
      if (testval!=BINFLAG) { 
	cerr << "Unrecognised binary matrix file format" << endl;
	return 2;
      }
    }

    // read matrix dimensions (num rows x num cols)
    unsigned int ival,nx,ny;
    fs.read((char*)&ival,sizeof(ival));
    // ignore the padding (reserved for future use)
    fs.read((char*)&ival,sizeof(ival));
    if (swapbytes) Swap_Nbytes(1,sizeof(ival),&ival);
    nx = ival;
    fs.read((char*)&ival,sizeof(ival));
    if (swapbytes) Swap_Nbytes(1,sizeof(ival),&ival);
    ny = ival;

    // set up and read matrix (rows fast, cols slow)
    double val;
    if ( (((unsigned int) mres.Ncols())<ny) || (((unsigned int) mres.Nrows())<nx) ) {
      mres.ReSize(nx,ny);
    }
    for (unsigned int y=1; y<=ny; y++) {
      for (unsigned int x=1; x<=nx; x++) {
	fs.read((char*)&val,sizeof(val));
	if (swapbytes) Swap_Nbytes(1,sizeof(val),&val);
	mres(x,y)=val;
      }
    }
    
    return 0;
  }


  // WRITE FUNCTIONS //


  int write_ascii_matrix(const string& filename, const Matrix& mat, 
			 int precision)
  {
    return write_ascii_matrix(mat, filename, precision);
  }

  int write_ascii_matrix(const Matrix& mat, const string& filename, 
			 int precision)
  {
    Tracer tr("write_ascii_matrix");
    if ( (filename.size()<1) ) return -1;
    ofstream fs(filename.c_str());
    if (!fs) { 
      cerr << "Could not open file " << filename << " for writing" << endl;
      return -1;
    }
    int retval = write_ascii_matrix(mat,fs,precision);
    fs.close();
    return retval;
  }

  int write_ascii_matrix(ofstream& fs, const Matrix& mat, 
			 int precision)
  {
    return write_ascii_matrix(mat, fs, precision);
  }
  
  int write_ascii_matrix(const Matrix& mat, ofstream& fs, int precision)
  {
    if (precision>0)  { 
      fs.setf(ios::scientific | ios::showpos); 
      fs.precision(precision); 
    } 
#ifdef PPC64	
    int n=0;
#endif
    for (int i=1; i<=mat.Nrows(); i++) {
      for (int j=1; j<=mat.Ncols(); j++) {
	fs << mat(i,j) << "  ";
#ifdef PPC64	
	if ((n++ % 50) == 0) fs.flush();
#endif
      }
      fs << endl;
    }
    return 0;
  }
  
  int write_vest(string p_fname, const Matrix& x, int precision)
     { return write_vest(x,p_fname,precision); } 

  int write_vest(const Matrix& x, string p_fname, int precision)
  {
    ofstream out;
    out.open(p_fname.c_str(), ios::out);
    
    if(!out)
      {
	cerr << "Unable to open " << p_fname << endl;
	return -1;
      }

    out << "! VEST-Waveform File" << endl;
    out << "/NumWaves\t" << x.Ncols() << endl;
    out << "/NumPoints\t" << x.Nrows() << endl;
    out << "/Skip" << endl;
    out << endl << "/Matrix" << endl;

    int retval = write_ascii_matrix(x, out, precision);

    out.close();

    return retval;
  }


  int write_binary_matrix(const Matrix& mat, const string& filename)
  {
    Tracer tr("write_binary_matrix");
    if ( (filename.size()<1) ) return -1;
    ofstream fs(filename.c_str(), ios::out | ios::binary);
    if (!fs) { 
      cerr << "Could not open file " << filename << " for writing" << endl;
      return -1;
    }
    int retval = write_binary_matrix(mat,fs);
    fs.close();
    return retval;
  }


  int write_binary_matrix(const Matrix& mat, ofstream& fs)
  {
    unsigned int ival, nx, ny;

    ival = BINFLAG;
    fs.write((char*)&ival,sizeof(ival));
    ival = 0;  // padding (reserved for future use)
    fs.write((char*)&ival,sizeof(ival));
    ival = mat.Nrows();
    fs.write((char*)&ival,sizeof(ival));
    ival = mat.Ncols();
    fs.write((char*)&ival,sizeof(ival));

    nx = mat.Nrows();
    ny = mat.Ncols();

    double val;
#ifdef PPC64	
    int n=0;
#endif
    for (unsigned int y=1; y<=ny; y++) {
      for (unsigned int x=1; x<=nx; x++) {
	val = mat(x,y);
	fs.write((char*)&val,sizeof(val));
#ifdef PPC64	
	if ((n++ % 50) == 0) fs.flush();
#endif
      }
    }

    return 0;
  }


  // General mathematical functions

  int round(int x) { return x; }

  int round(float x) 
    { 
      if (x>0.0) return ((int) (x+0.5));
      else       return ((int) (x-0.5));
    }

   int round(double x) 
   { 
     if (x>0.0) return ((int) (x+0.5));
     else       return ((int) (x-0.5));
   }  
  
  double rounddouble(double x){
    return ( floor(x+0.5));
  }
  int periodicclamp(int x, int x1, int x2)
   {
     if (x2<x1) return periodicclamp(x,x2,x1);
     int d = x2-x1+1;
     int xp = x-x1;
     if (xp>=0) {
       return (xp % d) + x1;
     } else {
       xp = xp + d + std::abs(xp/d)*d;
       assert(xp>0);
       return periodicclamp(xp + d + std::abs(xp/d)*d,x1,x2);
     }
   }

  ColumnVector cross(const ColumnVector& a, const ColumnVector& b)
    {
      Tracer tr("cross");
      ColumnVector ans(3);
      ans(1) = a(2)*b(3) - a(3)*b(2);
      ans(2) = a(3)*b(1) - a(1)*b(3);
      ans(3) = a(1)*b(2) - a(2)*b(1);
      return ans;
    }


  ColumnVector cross(const Real *a, const Real *b)
    {
      Tracer tr("cross");
      ColumnVector a1(3), b1(3);
      a1 << a;
      b1 << b;
      return cross(a1,b1);
    }


  double norm2(const ColumnVector& x)
    {
      return std::sqrt(x.SumSquare());
    }

  double norm2sq(double a, double b, double c)
    {
	return a*a + b*b + c*c;
    }

  float norm2sq(float a, float b, float c)
    {
	return a*a + b*b + c*c;
    }

  int diag(Matrix& m, const float diagvals[])
    {
      Tracer tr("diag");
      m=0.0;
      for (int j=1; j<=m.Nrows(); j++)
	m(j,j)=diagvals[j-1];
      return 0;
    }


  int diag(DiagonalMatrix& m, const ColumnVector& diagvals)
    {
      Tracer tr("diag");

      m.ReSize(diagvals.Nrows());
      m=0.0;
      for (int j=1; j<=diagvals.Nrows(); j++)
	m(j)=diagvals(j);
      return 0;
    }

  int diag(Matrix& m, const ColumnVector& diagvals)
    {
      Tracer tr("diag");
      
      m.ReSize(diagvals.Nrows(),diagvals.Nrows());
      m=0.0;
      for (int j=1; j<=diagvals.Nrows(); j++)
	m(j,j)=diagvals(j);
      return 0;
    }

  ReturnMatrix diag(const Matrix& Mat)
    {
      Tracer tr("diag");
      if(Mat.Ncols()==1){
	Matrix retmat(Mat.Nrows(),Mat.Nrows());
	diag(retmat,Mat);
	retmat.Release();
	return retmat;}
      else{
	int mindim = Min(Mat.Ncols(),Mat.Nrows());
	Matrix retmat(mindim,1);
	for(int ctr=1; ctr<=mindim;ctr++){
	  retmat(ctr,1)=Mat(ctr,ctr);
	}
	retmat.Release();
	return retmat;
      }
    }

  ReturnMatrix pinv(const Matrix& mat)
    {
      // calculates the psuedo-inverse using SVD
      // note that the right-pinv(x') = pinv(x).t()
      Tracer tr("pinv");
      DiagonalMatrix D;
      Matrix U, V;
      SVD(mat,D,U,V);
      float tol;
      tol = MaximumAbsoluteValue(D) * Max(mat.Nrows(),mat.Ncols()) * 1e-16;
      for (int n=1; n<=D.Nrows(); n++) {
	if (fabs(D(n,n))>tol) D(n,n) = 1.0/D(n,n);
	else D(n,n) = 0.0; // reduce the number of columns because too close to singular
      }
      Matrix pinv = V * D * U.t();
      pinv.Release();
      return pinv;
    }

  int rank(const Matrix& X)
    {
      // calculates the rank of matrix X
      Tracer tr("rank");

      DiagonalMatrix eigenvals;
      SVD(X,eigenvals);

      double tolerance = Max(X.Nrows(),X.Ncols()) * eigenvals.Maximum() * 1e-16;

      int therank=0;

      for(int i=0; i<eigenvals.Nrows(); i++)
	if (eigenvals(i+1)>tolerance)
	  therank++;

      // cout << "tolerance = " << tolerance << "\n" << "eigenvalues = " << eigenvals << "\n" << "rank = " << therank << endl;

      return therank;
    }


  ReturnMatrix sqrtaff(const Matrix& mat)
    {
      Tracer tr("sqrtaff");
      Matrix matnew(4,4), rot, id4;
      rot=IdentityMatrix(4);
      id4=IdentityMatrix(4);
      ColumnVector params(12), centre(3), trans(4);
      centre = 0.0;
      // Quaternion decomposition -> params(1..3) = sin(theta/2)*(unit_axis_vec)
      // Want a new quaternion : q = sin(theta/4)*(unit_axis_vec)
      // Therefore factor of conversion is: factor = sin(theta/4)/sin(theta/2)
      //      = 1/(2  * cos(theta/4))   which is calculated below
      //  NB: t = theta/2
      decompose_aff(params,mat,centre,rotmat2quat);
      double sint;
      sint = std::sqrt(params(1)*params(1) + params(2)*params(2) + 
		  params(3)*params(3));
      double t = asin(sint);
      double factor = 1.0/(2.0*cos(0.5*t));
      params(1) = factor * params(1);
      params(2) = factor * params(2);
      params(3) = factor * params(3);
      params(7) = std::sqrt(params(7));
      params(8) = std::sqrt(params(8));
      params(9) = std::sqrt(params(9));
      params(10) = 0.5*params(10);
      params(11) = 0.5*params(11);
      params(12) = 0.5*params(12);

      construct_rotmat_quat(params,3,rot,centre);
      rot(1,4) = 0.0;
      rot(2,4) = 0.0;
      rot(3,4) = 0.0;
      
      Matrix scale=IdentityMatrix(4);
      scale(1,1)=params(7);
      scale(2,2)=params(8);
      scale(3,3)=params(9);

      Matrix skew=IdentityMatrix(4);
      skew(1,2)=params(10);
      skew(1,3)=params(11);
      skew(2,3)=params(12);

      trans(1) = params(4);
      trans(2) = params(5);
      trans(3) = params(6);
      trans(4) = 1.0;

      // The translation, being independent of the 3x3 submatrix, is
      //  calculated so that it will be equal for each of the two 
      //  halves of the approximate square root 
      //  (i.e. matnew and mat*matnew.i() have exactly the same translation)
      ColumnVector th(4);
      th = (mat*scale.i()*skew.i()*rot.i() + id4).SubMatrix(1,3,1,3).i() 
	* trans.SubMatrix(1,3,1,1);

      matnew = rot*skew*scale;
      matnew(1,4) = th(1);
      matnew(2,4) = th(2);
      matnew(3,4) = th(3);

      matnew.Release();
      return matnew;
    }


  //------------------------------------------------------------------------//

  // Handy MATLAB-like functions
  
  void reshape(Matrix& r, const Matrix& m, int nrows, int ncols)
    {
      Tracer tr("reshape");
      if (nrows*ncols != m.Nrows() * m.Ncols() ) {
	cerr << "WARNING: cannot reshape " << m.Nrows() << "x"
	     << m.Ncols() << " matrix into " << nrows << "x" 
	     << ncols << endl;
	cerr << " Returning original matrix instead" << endl;
	r = m;
	return;
      }
      r.ReSize(nrows,ncols);
      int rr = 1, rc = 1;
      for (int mc=1; mc<=m.Ncols(); mc++) {
	for (int mr=1; mr<=m.Nrows(); mr++) {
	  r(rr,rc) = m(mr,mc);
	  rr++;
	  if (rr>nrows) {
	    rc++;
	    rr=1;
	  }
	}
      }
    }

  ReturnMatrix reshape(const Matrix& m, int nrows, int ncols)
    {
      Tracer tr("reshape");
     
      Matrix r; 
      
      reshape(r,m,nrows,ncols);

      r.Release();
      return r;
    }

  int addrow(Matrix& m, int ncols)
  {
    if (m.Nrows()==0) {
      Matrix mm(1,ncols);
      mm=0;
      m = mm;
    } else {
      Matrix mm(m.Nrows()+1,ncols);
      mm = 0;
      mm.SubMatrix(1,m.Nrows(),1,ncols) = m;
      m = mm;
    }
    return 0;
  }
  
  //------------------------------------------------------------------------//


  // Spatial transformation functions (rotations and affine transforms)


  int construct_rotmat_euler(const ColumnVector& params, int n, Matrix& aff,
			     const ColumnVector& centre)
    {
      Tracer tr("construct_rotmat_euler");
      ColumnVector angl(3);
      Matrix newaff(4,4);
      aff=IdentityMatrix(4);

      if (n<=0) return 0;
      // order of parameters is 3 rotation + 3 translation
      // angles are in radians 
      //  order of parameters is (Rx,Ry,Rz) and R = Rx.Ry.Rz
      angl=0.0;
      angl(1)=params(1);
      make_rot(angl,centre,newaff);
      aff = aff * newaff;
      if (n==1) return 0;

      angl=0.0;
      angl(2)=params(2);
      make_rot(angl,centre,newaff);
      aff = aff * newaff; 
      if (n==2) return 0;

      angl=0.0;
      angl(3)=params(3);
      make_rot(angl,centre,newaff);
      aff = aff * newaff;
      if (n==3) return 0;

      aff(1,4)+=params(4);
      if (n==4) return 0;
      aff(2,4)+=params(5);
      if (n==5) return 0;
      aff(3,4)+=params(6);
      if (n==6) return 0;

      return 1;
    }  

  int construct_rotmat_euler(const ColumnVector& params, int n, Matrix& aff)
    {
      Tracer tr("construct_rotmat_euler");
      ColumnVector centre(3);
      centre = 0.0;
      return construct_rotmat_euler(params,n,aff,centre);
    }


  int construct_rotmat_quat(const ColumnVector& params, int n, Matrix& aff,
			    const ColumnVector& centre)
    {
      Tracer tr("construct_rotmat_quat");
      aff=IdentityMatrix(4);

      if (n<=0) return 0;
      // order of parameters is 3 rotation (last 3 quaternion components) 
      //  + 3 translation

      if ((n>=1) && (n<3)) { cerr<<"Can only do 3 or more, not "<< n <<endl; }
      float w, w2 = 1.0 - Sqr(params(1)) - Sqr(params(2)) - Sqr(params(3));
      if (w2 < 0.0) {
	cerr << "Parameters do not form a valid axis - greater than unity\n";
	return -1;
      }
      w = std::sqrt(w2);
      float x=params(1), y=params(2), z=params(3);
      aff(1,1) = 1 - 2*y*y - 2*z*z;
      aff(2,2) = 1 - 2*x*x - 2*z*z;
      aff(3,3) = 1 - 2*x*x - 2*y*y;
      aff(1,2) = 2*x*y - 2*w*z;
      aff(2,1) = 2*x*y + 2*w*z;
      aff(1,3) = 2*x*z + 2*w*y;
      aff(3,1) = 2*x*z - 2*w*y;
      aff(2,3) = 2*y*z - 2*w*x;
      aff(3,2) = 2*y*z + 2*w*x;

      // Given Rotation matrix R:  x' = Rx + (I-R)*centre
      ColumnVector trans(3);
      trans = aff.SubMatrix(1,3,1,3)*centre;
      aff.SubMatrix(1,3,4,4) = centre - trans;

      aff(1,4) += params(4);
      if (n==4) return 0;
      aff(2,4) += params(5);
      if (n==5) return 0;
      aff(3,4) += params(6);
      if (n==6) return 0;

      return 1;
    }  

  int construct_rotmat_quat(const ColumnVector& params, int n, Matrix& aff)
    {
      Tracer tr("construct_rotmat_quat");
      ColumnVector centre(3);
      centre = 0.0;
      return construct_rotmat_quat(params,n,aff,centre);
    }

  int make_rot(const ColumnVector& angl, const ColumnVector& centre, 
	       Matrix& rot)
    {
      // Matrix rot must be 4x4; angl and orig must be length 3
      Tracer tr("make_rot");
      rot=IdentityMatrix(4);  // default return value
      float theta;
      theta = norm2(angl);
      if (theta<1e-8) {  // avoid round-off errors and return Identity
	return 0;
      }
      ColumnVector axis = angl/theta;
      ColumnVector x1(3), x2(3), x3(3);
      x1 = axis;
      x2(1) = -axis(2);  x2(2) = axis(1);  x2(3) = 0.0;
      if (norm2(x2)<=0.0) {
	x2(1) = 1.0;  x2(2) = 0.0;  x2(3) = 0.0;
      }
      x2 = x2/norm2(x2);
      x3 = cross(x1,x2);
      x3 = x3/norm2(x3);

      Matrix basischange(3,3);
      basischange.SubMatrix(1,3,1,1) = x2;
      basischange.SubMatrix(1,3,2,2) = x3;
      basischange.SubMatrix(1,3,3,3) = x1;

      Matrix rotcore=IdentityMatrix(3);
      rotcore(1,1)=cos(theta);
      rotcore(2,2)=cos(theta);
      rotcore(1,2)=sin(theta);
      rotcore(2,1)=-sin(theta);

      rot.SubMatrix(1,3,1,3) = basischange * rotcore * basischange.t();
  
      Matrix ident3=IdentityMatrix(3);
      ColumnVector trans(3);
      trans = (ident3 - rot.SubMatrix(1,3,1,3))*centre;
      rot.SubMatrix(1,3,4,4)=trans;
      return 0;
    }


  int getrotaxis(ColumnVector& axis, const Matrix& rotmat)
    {
      Tracer tr("getrotaxis");
      Matrix residuals(3,3);
      residuals = rotmat*rotmat.t() - IdentityMatrix(3);
      if (residuals.SumSquare() > 1e-4)
	{ cerr << "Failed orthogonality check!" << endl;  return -1; }
      Matrix u(3,3), v(3,3);
      DiagonalMatrix d(3);
      SVD(rotmat-IdentityMatrix(3),d,u,v);
      // return column of V corresponding to minimum value of |S|
      for (int i=1; i<=3; i++) {
	if (fabs(d(i))<1e-4)  axis = v.SubMatrix(1,3,i,i);
      }
      return 0;
    }  
  

  int rotmat2euler(ColumnVector& angles, const Matrix& rotmat)
    {
      // uses the convention that R = Rx.Ry.Rz
      Tracer tr("rotmat2euler");
      float cz, sz, cy, sy, cx, sx;
      cy = std::sqrt(Sqr(rotmat(1,1)) + Sqr(rotmat(1,2)));
      if (cy < 1e-4) {
	//cerr << "Cos y is too small - Gimbal lock condition..." << endl;
	cx = rotmat(2,2);
	sx = -rotmat(3,2);
	sy = -rotmat(1,3);
	angles(1) = atan2(sx,cx);
	angles(2) = atan2(sy,(float)0.0);
	angles(3) = 0.0;
      } else {
	// choose by convention that cy > 0
	// get the same rotation if: sy stays same & all other values swap sign
	cz = rotmat(1,1)/cy;
	sz = rotmat(1,2)/cy;
	cx = rotmat(3,3)/cy;
	sx = rotmat(2,3)/cy;
	sy = -rotmat(1,3);
	//atan2(sin,cos)  (defined as atan2(y,x))
	angles(1) = atan2(sx,cx);
	angles(2) = atan2(sy,cy);
	angles(3) = atan2(sz,cz);
      }
      return 0;
    }


  int rotmat2quat(ColumnVector& quaternion, const Matrix& rotmat)
    {
      Tracer tr("rotmat2quat");

      float trace = rotmat.SubMatrix(1,3,1,3).Trace();
	
      if (trace > 0) {
	float w = std::sqrt((trace + 1.0)/4.0);
	quaternion(1) = (rotmat(3,2) - rotmat(2,3))/(4.0*w);
	quaternion(2) = (rotmat(1,3) - rotmat(3,1))/(4.0*w);
	quaternion(3) = (rotmat(2,1) - rotmat(1,2))/(4.0*w);
      } else if ((rotmat(1,1) > rotmat(2,2)) && (rotmat(1,1) > rotmat(3,3))) {
	// first col case
	float s = std::sqrt(1.0 + rotmat(1,1) - rotmat(2,2) - rotmat(3,3)) * 2.0;
	quaternion(1) = 0.5 / s;
	quaternion(2) = (-rotmat(1,2) - rotmat(1,2)) / s;
	quaternion(3) = (-rotmat(1,3) - rotmat(3,1)) / s;
      } else if ((rotmat(2,2) > rotmat(1,1)) && (rotmat(2,2) > rotmat(3,3))) {
	// 2nd col case
	float s = std::sqrt(1.0 + rotmat(2,2) - rotmat(1,1) - rotmat(3,3)) * 2.0;
	quaternion(1) = (-rotmat(1,2) - rotmat(2,1)) / s;
	quaternion(2) = 0.5 / s;
	quaternion(3) = (-rotmat(2,3) - rotmat(3,2)) / s;
      } else if ((rotmat(3,3) > rotmat(1,1)) && (rotmat(3,3) > rotmat(2,2))) {
	// 3rd col case
	float s = std::sqrt(1.0 + rotmat(3,3) - rotmat(1,1) - rotmat(2,2)) * 2.0;
	quaternion(1) = (-rotmat(1,3) - rotmat(3,1)) / s;
	quaternion(2) = (-rotmat(2,3) - rotmat(3,2)) / s;
	quaternion(3) = 0.5 / s;
      }
      return 0;
    }


  int decompose_aff(ColumnVector& params, const Matrix& affmat, 
		    const ColumnVector& centre,
		    int (*rotmat2params)(ColumnVector& , const Matrix& ))
    {
      // decomposes using the convention: mat = rotmat * skew * scale
      // order of parameters is 3 rotation + 3 translation + 3 scales + 3 skews
      // angles are in radians
      Tracer tr("decompose_aff");
      if (params. Nrows() < 12)
	params.ReSize(12);
      if (rotmat2params==0)  
	{ 
	  cerr << "No rotmat2params function specified" << endl;  
	  return -1; 
	}
      ColumnVector x(3), y(3), z(3);
      Matrix aff3(3,3);
      aff3 = affmat.SubMatrix(1,3,1,3);
      x = affmat.SubMatrix(1,3,1,1);
      y = affmat.SubMatrix(1,3,2,2);
      z = affmat.SubMatrix(1,3,3,3);
      float sx, sy, sz, a, b, c;
      sx = norm2(x);
      sy = std::sqrt( dot(y,y) - (Sqr(dot(x,y)) / Sqr(sx)) );
      a = dot(x,y)/(sx*sy);
      ColumnVector x0(3), y0(3);
      x0 = x/sx;
      y0 = y/sy - a*x0;
      sz = std::sqrt(dot(z,z) - Sqr(dot(x0,z)) - Sqr(dot(y0,z)));
      b = dot(x0,z)/sz;
      c = dot(y0,z)/sz;
      params(7) = sx;  params(8) = sy;  params(9) = sz;
      Matrix scales(3,3);
      float diagvals[] = {sx,sy,sz};
      diag(scales,diagvals);
      Real skewvals[] = {1,a,b,0 , 0,1,c,0 , 0,0,1,0 , 0,0,0,1}; 
      Matrix skew(4,4);
      skew  << skewvals;
      params(10) = a;  params(11) = b;  params(12) = c;
      Matrix rotmat(3,3);
      rotmat = aff3 * scales.i() * (skew.SubMatrix(1,3,1,3)).i();
      ColumnVector transl(3);
      transl = affmat.SubMatrix(1,3,1,3)*centre + affmat.SubMatrix(1,3,4,4)
	         - centre;
      for (int i=1; i<=3; i++)  { params(i+3) = transl(i); }
      ColumnVector rotparams(3);
      (*rotmat2params)(rotparams,rotmat);
      for (int i=1; i<=3; i++)  { params(i) = rotparams(i); }
      return 0;
    }

  int decompose_aff(ColumnVector& params, const Matrix& affmat, 
		    int (*rotmat2params)(ColumnVector& , const Matrix& ))
    {
      Tracer tr("decompose_aff");
      ColumnVector centre(3);
      centre = 0.0;
      return decompose_aff(params,affmat,centre,rotmat2params);
    }



  int compose_aff(const ColumnVector& params, int n, const ColumnVector& centre,
		  Matrix& aff, 
		  int (*params2rotmat)(const ColumnVector& , int , Matrix& ,
			     const ColumnVector& ) )
    {
      Tracer tr("compose_aff");
      if (n<=0) return 0;
      // order of parameters is 3 rotation + 3 translation + 3 scales + 3 skews
      // angles are in radians

      (*params2rotmat)(params,n,aff,centre);
  
      if (n<=6)  return 0;
  
      Matrix scale=IdentityMatrix(4);
      if (n>=7) {
	scale(1,1)=params(7);
	if (n>=8) scale(2,2)=params(8);
	else      scale(2,2)=params(7);
	if (n>=9) scale(3,3)=params(9);
	else      scale(3,3)=params(7);
      }
      // fix the translation so that the centre is not moved
      ColumnVector strans(3);
      strans = centre - scale.SubMatrix(1,3,1,3)*centre;
      scale.SubMatrix(1,3,4,4) = strans;

      Matrix skew=IdentityMatrix(4);
      if (n>=10) {
	if (n>=10) skew(1,2)=params(10);
	if (n>=11) skew(1,3)=params(11);
	if (n>=12) skew(2,3)=params(12);
      }
      // fix the translation so that the centre is not moved
      ColumnVector ktrans(3);
      ktrans = centre - skew.SubMatrix(1,3,1,3)*centre;
      skew.SubMatrix(1,3,4,4) = ktrans;

      aff = aff * skew * scale;

      return 0;
    }


float rms_deviation(const Matrix& affmat1, const Matrix& affmat2, 
		    const ColumnVector& centre, const float rmax) 
{
  Tracer trcr("rms_deviation");
  Matrix isodiff(4,4), a1(4,4), a2(4,4);

  if ((affmat1.Nrows()==4) && (affmat1.Ncols()==4)) { a1=affmat1; }
  else if ((affmat1.Nrows()==3) && (affmat1.Ncols()==3)) { a1=IdentityMatrix(4); a1.SubMatrix(1,3,1,3)=affmat1; }
  else { cerr << "ERROR:: Can only calculate RMS deviation for 4x4 or 3x3 matrices" << endl; exit(-5); }

  if ((affmat2.Nrows()==4) && (affmat2.Ncols()==4)) { a2=affmat2; }
  else if ((affmat2.Nrows()==3) && (affmat2.Ncols()==3)) { a2=IdentityMatrix(4); a2.SubMatrix(1,3,1,3)=affmat2; }
  else { cerr << "ERROR:: Can only calculate RMS deviation for 4x4 or 3x3 matrices" << endl; exit(-5); }

  try {
    isodiff = a1*a2.i() - IdentityMatrix(4);
  } catch(...) {
    cerr << "RMS_DEVIATION ERROR:: Could not invert matrix" << endl;  
    exit(-5); 
  }
  Matrix adiff(3,3);
  adiff = isodiff.SubMatrix(1,3,1,3);
  ColumnVector tr(3);
  tr = isodiff.SubMatrix(1,3,4,4) + adiff*centre;
  float rms = std::sqrt( (tr.t() * tr).AsScalar() + 
		    (rmax*rmax/5.0)*Trace(adiff.t()*adiff) );
  return rms;
}


float rms_deviation(const Matrix& affmat1, const Matrix& affmat2, 
		    const float rmax) 
{
  ColumnVector centre(3);
  centre = 0;
  return rms_deviation(affmat1,affmat2,centre,rmax);
}


  // helper function - calls nifti, but with FSL default case

Matrix Mat44ToNewmat(mat44 m)
{
  Matrix r(4,4);

  for(unsigned short i = 0; i < 4; ++i)
    for(unsigned short j = 0; j < 4; ++j)
      r(i+1, j+1) = m.m[i][j];
      
  return r;
}

mat44 NewmatToMat44(const Matrix& m)
{
  mat44 r;

  for(unsigned short i = 0; i < 4; ++i)
    for(unsigned short j = 0; j < 4; ++j)
      r.m[i][j] = m(i+1, j+1);

  return r;
}

void get_axis_orientations(const Matrix& sform_mat, int sform_code,
			   const Matrix& qform_mat, int qform_code,
			   int& icode, int& jcode, int& kcode)
{
  Matrix vox2mm(4,4);
  if (sform_code!=NIFTI_XFORM_UNKNOWN) {
    vox2mm = sform_mat;
  } else if (qform_code!=NIFTI_XFORM_UNKNOWN) {
    vox2mm = qform_mat;
  } else {
    // ideally should be sampling_mat(), but for orientation it doesn't matter
    vox2mm = IdentityMatrix(4);
    vox2mm(1,1) = -vox2mm(1,1);
  }
  mat44 v2mm;
  for (int ii=0; ii<4; ii++) { for (int jj=0; jj<4; jj++) {
      v2mm.m[ii][jj] = vox2mm(ii+1,jj+1);
    } }
  nifti_mat44_to_orientation(v2mm,&icode,&jcode,&kcode);
}

 
Matrix mat44_to_newmat(mat44 inmat)
{
  Matrix retmat(4,4);
  for (int ii=0; ii<4; ii++) {
    for (int jj=0; jj<4; jj++) {
      retmat(ii+1,jj+1) = inmat.m[ii][jj];
    }
  }
  return retmat;
}
 	 
mat44 newmat_to_mat44(const Matrix& inmat)
{
  mat44 retmat;
  for (int ii=0; ii<4; ii++) {
    for (int jj=0; jj<4; jj++) {
      retmat.m[ii][jj] = inmat(ii+1,jj+1);
    }
  }
  return retmat;
}
 	 
// Matlab style functions for percentiles, quantiles and median
// AUG 06 CB

ColumnVector seq(const int size)
{
  ColumnVector outputVector(size);
  for(int i=1; i<=size; i++)
    outputVector(i) = i;
  return outputVector;
}

float interp1(const ColumnVector& x, const ColumnVector& y, float xi) 
// Look-up function for data table defined by x, y
// Returns the values yi at xi using linear interpolation
// Assumes that x is sorted in ascending order
{
  
  float ans;
  if(xi >= x.Maximum()) 
    ans = y(x.Nrows());
  else
    if(xi <= x.Minimum()) 
      ans = y(1); 
    else{
      int ind=1;
      while(xi >= x(ind))
	ind++;      
      float xa = x(ind-1), xb = x(ind), ya = y(ind-1), yb = y(ind);
      ans = ya + (xi - xa)/(xb - xa) * (yb - ya);
    }
  return ans;
}


float quantile(const ColumnVector& in, int which)
{
  float p;
  switch (which)
    {  
    case 0 : p =  0.0; break;
    case 1 : p = 25.0; break;
    case 2 : p = 50.0; break; 
    case 3 : p = 75.0; break;
    case 4 : p =100.0; break;
    default: p =  0.0; 
    }

  return percentile(in,p);
}

float percentile(const ColumnVector& in, float p)
{
  ColumnVector y = in;
  SortAscending(y);
  int num = y.Nrows();

  ColumnVector xx,yy,sequence,a(1),b(1),c(1),d(1);
  sequence = 100*(seq(num)-0.5)/num; a << y(1); b << y(num); c = 0; d = 100;
  xx = (c & sequence & d);
  yy = (a & y & b);
  
  return interp1(xx,yy,p);
}

ReturnMatrix quantile(const Matrix& in, int which)
{
  int num = in.Ncols();
  Matrix res(1,num);
  for (int ctr=1; ctr<=num; ctr++){
    ColumnVector tmp = in.Column(ctr);
    res(1,ctr) = quantile(tmp,which);
  }
  res.Release();
  return res;
}

ReturnMatrix  percentile(const Matrix& in, float p)
{
  int num = in.Ncols();
  Matrix res(1,num);
  for (int ctr=1; ctr<=num; ctr++){
    ColumnVector tmp = in.Column(ctr);
    res(1,ctr) = percentile(tmp,p);
  }
  res.Release();
  return res;
}



void cart2sph(const ColumnVector& dir, float& th, float& ph)
{
  float mag=sqrt(dir(1)*dir(1)+dir(2)*dir(2)+dir(3)*dir(3));
  if(mag==0){
    ph=M_PI/2;
    th=M_PI/2;
  }
  else{

    if(dir(1)==0 && dir(2)>=0) ph=M_PI/2;
    else if(dir(1)==0 && dir(2)<0) ph=-M_PI/2;
    else if(dir(1)>0) ph=atan(dir(2)/dir(1));
    else if(dir(2)>0) ph=atan(dir(2)/dir(1))+M_PI;
    else ph=atan(dir(2)/dir(1))-M_PI;
    
    if(dir(3)==0) th=M_PI/2;
    else if(dir(3)>0) th=atan(sqrt(dir(1)*dir(1)+dir(2)*dir(2))/dir(3));
    else th=atan(sqrt(dir(1)*dir(1)+dir(2)*dir(2))/dir(3))+M_PI;
  }
}



void cart2sph(const Matrix& dir,ColumnVector& th,ColumnVector& ph)
{
  if(th.Nrows()!=dir.Ncols()){
    th.ReSize(dir.Ncols());
  }

  if(ph.Nrows()!=dir.Ncols()){
    ph.ReSize(dir.Ncols());
  }

  for (int i=1;i<=dir.Ncols();i++) {
    float mag=sqrt(dir(1,i)*dir(1,i)+dir(2,i)*dir(2,i)+dir(3,i)*dir(3,i));
    if(mag==0){
      ph(i)=M_PI/2;
      th(i)=M_PI/2;
    }
    else{
      if(dir(1,i)==0 && dir(2,i)>=0) ph(i)=M_PI/2;
      else if(dir(1,i)==0 && dir(2,i)<0) ph(i)=-M_PI/2;
      else if(dir(1,i)>0) ph(i)=atan(dir(2,i)/dir(1,i));
      else if(dir(2,i)>0) ph(i)=atan(dir(2,i)/dir(1,i))+M_PI;
      else ph(i)=atan(dir(2,i)/dir(1,i))-M_PI;

      if(dir(3,i)==0) th(i)=M_PI/2;
      else if(dir(3,i)>0) th(i)=atan(sqrt(dir(1,i)*dir(1,i)+dir(2,i)*dir(2,i))/dir(3,i));
      else th(i)=atan(sqrt(dir(1,i)*dir(1,i)+dir(2,i)*dir(2,i))/dir(3,i))+M_PI;

    }
  }
}

// added by SJ
void cart2sph(const vector<ColumnVector>& dir,ColumnVector& th,ColumnVector& ph)
{
  if(th.Nrows()!=(int)dir.size()){
    th.ReSize(dir.size());
  }

  if(ph.Nrows()!=(int)dir.size()){
    ph.ReSize(dir.size());
  }
  //double _2pi=2*M_PI;
  double _pi2=M_PI/2;

  int j=1;
  for (unsigned int i=0;i<dir.size();i++) {
	float mag=std::sqrt(dir[i](1)*dir[i](1)+dir[i](2)*dir[i](2)+dir[i](3)*dir[i](3));
    if(mag==0){
      ph(j)=_pi2;
      th(j)=_pi2;
    }
    else{
      if(dir[i](1)==0 && dir[i](2)>=0) ph(j)=_pi2;
      else if(dir[i](1)==0 && dir[i](2)<0) ph(j)=-_pi2;
      else if(dir[i](1)>0) ph(j)=std::atan(dir[i](2)/dir[i](1));
      else if(dir[i](2)>0) ph(j)=std::atan(dir[i](2)/dir[i](1))+M_PI;
      else ph(j)=std::atan(dir[i](2)/dir[i](1))-M_PI;

      //ph(j)=fmod(ph(j),_2pi);

      if(dir[i](3)==0) th(j)=_pi2;
      else if(dir[i](3)>0) th(j)=std::atan(std::sqrt(dir[i](1)*dir[i](1)+dir[i](2)*dir[i](2))/dir[i](3));
      else th(j)=std::atan(std::sqrt(dir[i](1)*dir[i](1)+dir[i](2)*dir[i](2))/dir[i](3))+M_PI;
      
      //th(j)=fmod(th(j),M_PI);

    }
    j++;
  }
}

// Added by CFB   --- Matlab style Matrix functions
 
ReturnMatrix ones(const int dim1, const int dim2)
{ 
  int tdim = dim2;
  if(tdim<0){tdim=dim1;}
  Matrix res(dim1,tdim); res = 1.0;
  res.Release();
  return res;
}

ReturnMatrix zeros(const int dim1, const int dim2)
{ 
  int tdim = dim2;
  if(tdim<0){tdim=dim1;}
  Matrix res(dim1,tdim); res = 0.0;
  res.Release();
  return res;
}

ReturnMatrix repmat(const Matrix &mat, const int rows, const int cols)
{
  Matrix res = mat;
  for(int ctr = 1; ctr < cols; ctr++){res |= mat;}
  Matrix tmpres = res;
  for(int ctr = 1; ctr < rows; ctr++){res &= tmpres;}
  res.Release();
  return res;
}

ReturnMatrix dist2(const Matrix &mat1, const Matrix &mat2)
{
  Matrix res(mat1.Ncols(),mat2.Ncols());
  for(int ctr1 = 1; ctr1 <= mat1.Ncols(); ctr1++)
    for(int ctr2 =1; ctr2 <= mat2.Ncols(); ctr2++)
      {
	ColumnVector tmp;
	tmp=mat1.Column(ctr1)-mat2.Column(ctr2);
	res(ctr1,ctr2) = std::sqrt(tmp.SumSquare());
      }
  res.Release();
  return res;
}

ReturnMatrix abs(const Matrix& mat)
{
  Matrix res = mat;
  for (int mc=1; mc<=mat.Ncols(); mc++) {
    for (int mr=1; mr<=mat.Nrows(); mr++) {
      res(mr,mc)=std::abs(res(mr,mc));
    }
  }
  res.Release();
  return res;
}

ReturnMatrix sqrt(const Matrix& mat)
{
  Matrix res = mat;
  bool neg_flag = false;
  for (int mc=1; mc<=mat.Ncols(); mc++) {
    for (int mr=1; mr<=mat.Nrows(); mr++) {
      if(res(mr,mc)<0){ neg_flag = true; }
      res(mr,mc)=std::sqrt(std::abs(res(mr,mc)));
    }
  }
  if(neg_flag){
    //cerr << " Matrix contained negative elements" << endl;
    //cerr << " return sqrt(abs(X)) instead" << endl;
  }
  res.Release();
  return res;
}

ReturnMatrix sqrtm(const Matrix& mat)
{
	Matrix res, tmpU, tmpV;
	DiagonalMatrix tmpD;
	SVD(mat, tmpD, tmpU, tmpV);
	res = tmpU*sqrt(tmpD)*tmpV.t();
	res.Release();
	return res;
}

ReturnMatrix log(const Matrix& mat)
{
  Matrix res = mat;
  bool neg_flag = false;
  for (int mc=1; mc<=mat.Ncols(); mc++) {
    for (int mr=1; mr<=mat.Nrows(); mr++) {
      if(res(mr,mc)<0){ neg_flag = true; }
      res(mr,mc)=std::log(std::abs(res(mr,mc)));
    }
  }
  if(neg_flag){
    //  cerr << " Matrix contained negative elements" << endl;
    //  cerr << " return log(abs(X)) instead" << endl;
  }
  res.Release();
  return res; 
}

ReturnMatrix exp(const Matrix& mat)
{
  Matrix res = mat;
  for (int mc=1; mc<=mat.Ncols(); mc++) {
    for (int mr=1; mr<=mat.Nrows(); mr++) {
      res(mr,mc)=std::exp(res(mr,mc));
    }
  }
  res.Release();
  return res;
}

  // optimised code for calculating matrix exponential
ReturnMatrix expm(const Matrix& mat){
  float nmat = sum(mat).Maximum();
  int nc=mat.Ncols(),nr=mat.Nrows();
  Matrix res(nr,nc);
  IdentityMatrix id(nr);
  Matrix U(nr,nc),V(nr,nc);

  if(nmat <= 1.495585217958292e-002){ // m=3
    Matrix mat2(nr,nc);
    mat2=mat*mat;
    U = mat*(mat2+60.0*id);
    V = 12.0*mat2+120.0*id;
    res = (-U+V).i()*(U+V);
  }
  else if(nmat <= 2.539398330063230e-001){ // m=5
    Matrix mat2(nr,nc),mat4(nr,nc);
    mat2=mat*mat;mat4=mat2*mat2;
    U = mat*(mat4+420.0*mat2+15120.0*id);
    V = 30.0*mat4+3360.0*mat2+30240.0*id;
    res = (-U+V).i()*(U+V);
  }
  else if(nmat <= 9.504178996162932e-001){ // m=7
    Matrix mat2(nr,nc),mat4(nr,nc),mat6(nr,nc);
    mat2=mat*mat;mat4=mat2*mat2,mat6=mat4*mat2;
    U = mat*(mat6+1512.0*mat4+277200.0*mat2+8648640.0*id);
    V = 56.0*mat6+25200.0*mat4+1995840.0*mat2+17297280.0*id;
    res = (-U+V).i()*(U+V);
  }
  else if(nmat <= 2.097847961257068e+000){
    Matrix mat2(nr,nc),mat4(nr,nc),mat6(nr,nc),mat8(nr,nc);
    mat2=mat*mat;mat4=mat2*mat2,mat6=mat4*mat2,mat8=mat6*mat2;
    U = mat*(mat8+3960.0*mat6+2162160.0*mat4+302702400.0*mat2+8821612800.0*id);
    V = 90.0*mat8+110880.0*mat6+30270240.0*mat4+2075673600.0*mat2+17643225600.0*id;
    res = (-U+V).i()*(U+V);
  }
  else if(nmat <= 5.371920351148152e+000){
    Matrix mat2(nr,nc),mat4(nr,nc),mat6(nr,nc);
    mat2=mat*mat;mat4=mat2*mat2,mat6=mat4*mat2;
    U = mat*(mat6*(mat6+16380.0*mat4+40840800.0*mat2)+
	     +33522128640.0*mat6+10559470521600.0*mat4+1187353796428800.0*mat2+32382376266240000.0*id);
    V = mat6*(182.0*mat6+960960.0*mat4+1323241920.0*mat2)
      + 670442572800.0*mat6+129060195264000.0*mat4+7771770303897600.0*mat2+64764752532480000.0*id;
    res = (-U+V).i()*(U+V);
  }
  else{
    double t;int s;
    t = frexp(nmat/5.371920351148152,&s);
    if(t==0.5) s--;
    t = std::pow(2.0,s);
    res = (mat/t);
    Matrix mat2(nr,nc),mat4(nr,nc),mat6(nr,nc);
    mat2=res*res;mat4=mat2*mat2,mat6=mat4*mat2;
    U = res*(mat6*(mat6+16380*mat4+40840800*mat2)+
	     +33522128640.0*mat6+10559470521600.0*mat4+1187353796428800.0*mat2+32382376266240000.0*id);
    V = mat6*(182.0*mat6+960960.0*mat4+1323241920.0*mat2)
      + 670442572800.0*mat6+129060195264000.0*mat4+7771770303897600.0*mat2+64764752532480000.0*id;
    res = (-U+V).i()*(U+V);
    for(int i=1;i<=s;i++)
      res = res*res;
  }

  res.Release();
  return res;
}


ReturnMatrix tanh(const Matrix& mat)
{
  Matrix res = mat;
  for (int mc=1; mc<=mat.Ncols(); mc++) {
    for (int mr=1; mr<=mat.Nrows(); mr++) {
      res(mr,mc)=std::tanh(res(mr,mc));
    }
  }
  res.Release();
  return res;
}

ReturnMatrix pow(const Matrix& mat, const double exp)
{
  Matrix res = mat;
  for (int mc=1; mc<=mat.Ncols(); mc++) {
    for (int mr=1; mr<=mat.Nrows(); mr++) {
      res(mr,mc)=std::pow(res(mr,mc),exp);
    }
  }
  res.Release();
  return res;
}

ReturnMatrix max(const Matrix& mat)
{
  Matrix res;
  if(mat.Nrows()>1){
    res=zeros(1,mat.Ncols());
    res=mat.Row(1);
    for(int mc=1; mc<=mat.Ncols();mc++){
      for(int mr=2; mr<=mat.Nrows();mr++){
	if(mat(mr,mc)>res(1,mc)){res(1,mc)=mat(mr,mc);}
      }
    }
  }
  else{
    res=zeros(1);
    res=mat(1,1);
    for(int mc=2; mc<=mat.Ncols(); mc++){
      if(mat(1,mc)>res(1,1)){res(1,1)=mat(1,mc);}
    }
  }
  res.Release();
  return res;
}

ReturnMatrix max(const Matrix& mat,ColumnVector& index)
{
  index.ReSize(mat.Nrows());
  index=1;
  Matrix res;
  if(mat.Nrows()>1){
    res=zeros(1,mat.Ncols());
    res=mat.Row(1);
    for(int mc=1; mc<=mat.Ncols();mc++){
      for(int mr=2; mr<=mat.Nrows();mr++){
	if(mat(mr,mc)>res(1,mc))
	  {
	    res(1,mc)=mat(mr,mc);
	    index(mr)=mc;
	  }
      }
    }
  }
  else{
    res=zeros(1);
    res=mat(1,1);
    for(int mc=2; mc<=mat.Ncols(); mc++){
      if(mat(1,mc)>res(1,1))
	{
	  res(1,1)=mat(1,mc);
	  index(1)=mc;
	}
    }
  }
  res.Release();
  return res;
}

ReturnMatrix min(const Matrix& mat)
{
  Matrix res;
  if(mat.Nrows()>1){
    res=zeros(1,mat.Ncols());
    res=mat.Row(1);
    for(int mc=1; mc<=mat.Ncols();mc++){
      for(int mr=2; mr<=mat.Nrows();mr++){
	if(mat(mr,mc)<res(1,mc)){res(1,mc)=mat(mr,mc);}
      }
    }
  }
  else{
    res=zeros(1);
    res=mat(1,1);
    for(int mc=2; mc<=mat.Ncols(); mc++){
      if(mat(1,mc)<res(1,1)){res(1,1)=mat(1,mc);}
    }
  }
  res.Release();
  return res;
}
  

ReturnMatrix sum(const Matrix& mat, const int dim)
{
  Matrix tmp;

  if (dim == 1) {tmp=mat;}
  else {tmp=mat.t();}
  Matrix res(1,tmp.Ncols());
  res = 0.0;  
  for (int mc=1; mc<=tmp.Ncols(); mc++) {
    for (int mr=1; mr<=tmp.Nrows(); mr++) {
      res(1,mc) += tmp(mr,mc);
    }
  }
  if (!(dim == 1)) {res=res.t();}
  res.Release();
  return res;
}

ReturnMatrix mean(const Matrix& mat, const int dim)
{
  Matrix tmp;
  if (dim == 1) {tmp=mat;}
  else {tmp=mat.t();}

  int N = tmp.Nrows();

  Matrix res(1,tmp.Ncols());
  res = 0.0;  
  for (int mc=1; mc<=tmp.Ncols(); mc++) {
    for (int mr=1; mr<=tmp.Nrows(); mr++) {
      res(1,mc) += tmp(mr,mc)/N;
    }
  }
  if (!(dim == 1)) {res=res.t();}

  res.Release();
  return res;
}

ReturnMatrix var(const Matrix& mat, const int dim)
{
  Matrix tmp;
  if (dim == 1) {tmp=mat;}
  else {tmp=mat.t();}
  int N = tmp.Nrows();
  Matrix res(1,tmp.Ncols());
  res = 0.0;

  if(N>1){    
    tmp -= ones(tmp.Nrows(),1)*mean(tmp,1);   
    for (int mc=1; mc<=tmp.Ncols(); mc++) 
      for (int mr=1; mr<=tmp.Nrows(); mr++) 
        res(1,mc) += tmp(mr,mc) / (N-1) * tmp(mr,mc);
  }

  if (!(dim == 1)) {res=res.t();}
  res.Release();
  return res;
}

ReturnMatrix stdev(const Matrix& mat, const int dim)
{
  return sqrt(var(mat,dim));
}

ReturnMatrix gt(const Matrix& mat1,const Matrix& mat2)
{
  int ctrcol = std::min(mat1.Ncols(),mat2.Ncols());
  int ctrrow = std::min(mat1.Nrows(),mat2.Nrows());
  Matrix res(ctrrow,ctrcol);
  res=0.0;

  for (int ctr1 = 1; ctr1 <= ctrrow; ctr1++) {
    for (int ctr2 =1; ctr2 <= ctrcol; ctr2++) {
      if( mat1(ctr1,ctr2) > mat2(ctr1,ctr2)){
	res(ctr1,ctr2) = 1.0;
      }
    }
  }

  res.Release();
  return res;
}


ReturnMatrix lt(const Matrix& mat1,const Matrix& mat2)
{
  int ctrcol = std::min(mat1.Ncols(),mat2.Ncols());
  int ctrrow = std::min(mat1.Nrows(),mat2.Nrows());
  Matrix res(ctrrow,ctrcol);
  res=0.0;

  for (int ctr1 = 1; ctr1 <= ctrrow; ctr1++) {
    for (int ctr2 =1; ctr2 <= ctrcol; ctr2++) {
      if( mat1(ctr1,ctr2) < mat2(ctr1,ctr2)){
	res(ctr1,ctr2) = 1.0;
      }
    }
  }

  res.Release();
  return res;
}


ReturnMatrix geqt(const Matrix& mat1,const Matrix& mat2) 
{
  int ctrcol = std::min(mat1.Ncols(),mat2.Ncols());
  int ctrrow = std::min(mat1.Nrows(),mat2.Nrows());
  Matrix res(ctrrow,ctrcol);
  res=0.0;

  for (int ctr1 = 1; ctr1 <= ctrrow; ctr1++) {
    for (int ctr2 =1; ctr2 <= ctrcol; ctr2++) {
      if( mat1(ctr1,ctr2) >= mat2(ctr1,ctr2)){
	res(ctr1,ctr2) = 1.0;
      }
    }
  }

  res.Release();
  return res;
}
ReturnMatrix geqt(const Matrix& mat,const float a) 
{
  int ncols = mat.Ncols();
  int nrows = mat.Nrows();
  Matrix res(nrows,ncols);
  res=0.0;

  for (int ctr1 = 1; ctr1 <= nrows; ctr1++) {
    for (int ctr2 =1; ctr2 <= ncols; ctr2++) {
      if( mat(ctr1,ctr2) >= a){
	res(ctr1,ctr2) = 1.0;
      }
    }
  }

  res.Release();
  return res;
}

ReturnMatrix leqt(const Matrix& mat1,const Matrix& mat2) 
{
  int ctrcol = std::min(mat1.Ncols(),mat2.Ncols());
  int ctrrow = std::min(mat1.Nrows(),mat2.Nrows());
  Matrix res(ctrrow,ctrcol);
  res=0.0;

  for (int ctr1 = 1; ctr1 <= ctrrow; ctr1++) {
    for (int ctr2 =1; ctr2 <= ctrcol; ctr2++) {
      if( mat1(ctr1,ctr2) <= mat2(ctr1,ctr2)){
	res(ctr1,ctr2) = 1.0;
      }
    }
  }

  res.Release();
  return res;
}


ReturnMatrix eq(const Matrix& mat1,const Matrix& mat2)
{
  int ctrcol = std::min(mat1.Ncols(),mat2.Ncols());
  int ctrrow = std::min(mat1.Nrows(),mat2.Nrows());
  Matrix res(ctrrow,ctrcol);
  res=0.0;

  for (int ctr1 = 1; ctr1 <= ctrrow; ctr1++) {
    for (int ctr2 =1; ctr2 <= ctrcol; ctr2++) {
      if( mat1(ctr1,ctr2) == mat2(ctr1,ctr2)){
	res(ctr1,ctr2) = 1.0;
      }
    }
  }

  res.Release();
  return res;
}


ReturnMatrix neq(const Matrix& mat1,const Matrix& mat2) 
{
  int ctrcol = std::min(mat1.Ncols(),mat2.Ncols());
  int ctrrow = std::min(mat1.Nrows(),mat2.Nrows());
  Matrix res(ctrrow,ctrcol);
  res=0.0;

  for (int ctr1 = 1; ctr1 <= ctrrow; ctr1++) {
    for (int ctr2 =1; ctr2 <= ctrcol; ctr2++) {
      if( mat1(ctr1,ctr2) != mat2(ctr1,ctr2)){
	res(ctr1,ctr2) = 1.0;
      }
    }
  }

  res.Release();
  return res;
}

ReturnMatrix SD(const Matrix& mat1,const Matrix& mat2) 
{
  if((mat1.Nrows() != mat2.Nrows()) ||
     (mat1.Ncols() != mat2.Ncols()) ){
    cerr <<"MISCMATHS::SD - matrices are of different dimensions"<<endl;
    exit(-1);
  }
  Matrix ret(mat1.Nrows(),mat1.Ncols());
  for (int r = 1; r <= mat1.Nrows(); r++) {
    for (int c =1; c <= mat1.Ncols(); c++) {
      if( mat2(r,c)==0)
	ret(r,c)=0;
      else
	ret(r,c) = mat1(r,c)/mat2(r,c);
    }
  }

  ret.Release();
  return ret;
}

//Deprecate?
ReturnMatrix vox_to_vox(const ColumnVector& xyz1,const ColumnVector& dims1,const ColumnVector& dims2,const Matrix& xfm){
  ColumnVector xyz1_mm(4),xyz2_mm,xyz2(3);
  xyz1_mm<<xyz1(1)*dims1(1)<<xyz1(2)*dims1(2)<<xyz1(3)*dims1(3)<<1;
  xyz2_mm=xfm*xyz1_mm;
  xyz2_mm=xyz2_mm/xyz2_mm(4);
  xyz2<<xyz2_mm(1)/dims2(1)<<xyz2_mm(2)/dims2(2)<<xyz2_mm(3)/dims2(3);
  xyz2.Release();
  return xyz2;
}

//Deprecate?
ReturnMatrix mni_to_imgvox(const ColumnVector& mni,const ColumnVector& mni_origin,const Matrix& mni2img, const ColumnVector& img_dims){
  ColumnVector mni_new_origin(4),img_mm;//homogeneous
  ColumnVector img_vox(3);
  mni_new_origin<<mni(1)+mni_origin(1)<<mni(2)+mni_origin(2)<<mni(3)+mni_origin(3)<<1;
  img_mm=mni2img*mni_new_origin;
  img_vox<<img_mm(1)/img_dims(1)<<img_mm(2)/img_dims(2)<<img_mm(3)/img_dims(3);
  img_vox.Release();
  return img_vox;
}


ReturnMatrix remmean(const Matrix& mat, const int dim)
{ 
  Matrix res;
  if (dim == 1) {res=mat;}
  else {res=mat.t();}

  Matrix Mean;
  Mean = mean(res);

  for (int ctr = 1; ctr <= res.Nrows(); ctr++) {
    for (int ctr2 =1; ctr2 <= res.Ncols(); ctr2++) {
      res(ctr,ctr2)-=Mean(1,ctr2);
    }
  }
  if (dim != 1) {res=res.t();}
  res.Release();
  return res;
}


void remmean(const Matrix& mat, Matrix& demeanedmat, Matrix& Mean,  const int dim)
{ 
  if (dim == 1) {demeanedmat=mat;}
  else {demeanedmat=mat.t();}

  Mean = mean(demeanedmat);

  for (int ctr = 1; ctr <= demeanedmat.Nrows(); ctr++) {
    for (int ctr2 =1; ctr2 <= demeanedmat.Ncols(); ctr2++) {
      demeanedmat(ctr,ctr2)-=Mean(1,ctr2);
    }
  }
  if (dim != 1){demeanedmat = demeanedmat.t();Mean = Mean.t();}
}

ReturnMatrix cov(const Matrix& mat, const int norm)
{ 
  SymmetricMatrix res;
  Matrix tmp;
  int N;
  tmp=remmean(mat);
  if (norm == 1) {N = mat.Nrows();}
  else {N = mat.Nrows()-1;}  
  res << tmp.t()*tmp;
  res = res/N;

  res.Release();
  return res; 
}

ReturnMatrix corrcoef(const Matrix& mat, const int norm)
{ 
  SymmetricMatrix res;
  SymmetricMatrix C;
  C = cov(mat,norm);
  Matrix D;
  D=diag(C);
  D=pow(sqrt(D*D.t()),-1);
  res << SP(C,D);
  res.Release();
  return res; 
}

ReturnMatrix flipud(const Matrix& mat)
{
  Matrix rmat(mat.Nrows(),mat.Ncols());
  for(int j=1;j<=mat.Ncols();j++)
    for(int i=1;i<=mat.Nrows();i++)
      rmat(i,j)=mat(mat.Nrows()-i+1,j);
  rmat.Release();
  return rmat;
}

ReturnMatrix fliplr(const Matrix& mat)
{
  Matrix rmat(mat.Nrows(),mat.Ncols());
  for(int j=1;j<=mat.Ncols();j++)
    for(int i=1;i<=mat.Nrows();i++)
      rmat(i,j)=mat(i,mat.Ncols()-j+1);
  rmat.Release();
  return rmat;
}


void symm_orth(Matrix &Mat)
{
  SymmetricMatrix Metric;
  Metric << Mat.t()*Mat;
  Metric = Metric.i();
  Matrix tmpE;
  DiagonalMatrix tmpD;
  EigenValues(Metric,tmpD,tmpE);
  Mat = Mat * tmpE * sqrt(abs(tmpD)) * tmpE.t();
}

void powerspectrum(const Matrix &Mat1, Matrix &Result, bool useLog)
  //calculates the powerspectrum for every column of Mat1
{
  Matrix res;
  for(int ctr=1; ctr <= Mat1.Ncols(); ctr++)
    {
      ColumnVector tmpCol;
      tmpCol=Mat1.Column(ctr);
      ColumnVector FtmpCol_real;
      ColumnVector FtmpCol_imag;
      ColumnVector tmpPow;
      
      RealFFT(tmpCol,FtmpCol_real,FtmpCol_imag);
      tmpPow = pow(FtmpCol_real,2)+pow(FtmpCol_imag,2);
      tmpPow = tmpPow.Rows(2,tmpPow.Nrows());
      if(useLog){tmpPow = log(tmpPow);}
      if(res.Storage()==0){res= tmpPow;}else{res|=tmpPow;}
    }
    Result=res;
}


void element_mod_n(Matrix& Mat,double n)
{
 //represent each element in modulo n (useful for wrapping phases (n=2*M_PI))

  double tmp;
  for( int j=1;j<=Mat.Ncols();j++){
    for( int i=1;i<=Mat.Nrows();i++){

      while( !( (Mat(i,j)>0) && (Mat(i,j)<n) ) ){ 
	tmp = ( Mat(i,j) - rounddouble(Mat(i,j)/n)*n );
	Mat(i,j)= tmp > 0 ? tmp : tmp + n;
      }    
      
    }
    
  }
  
}

int nextpow2(int n)
{
  return (int)pow(2,ceil(log(float(n))/log(float(2))));
}

void xcorr(const Matrix& p_ts, Matrix& ret, int lag, int p_zeropad)
{
  Tracer tr("MISCMATHS::xcorr");

  int sizeTS = p_ts.Nrows();
  int numTS = p_ts.Ncols();
  
  if(p_zeropad == 0)
    p_zeropad = sizeTS;
  if(lag == 0)
    lag = sizeTS;

  ColumnVector x(p_zeropad);
  x = 0;
  ColumnVector fft_real;
  ColumnVector fft_im;
  ColumnVector dummy(p_zeropad);
  ColumnVector dummy2;
  dummy = 0;
  ColumnVector realifft(p_zeropad);
  ret.ReSize(lag,numTS);
  ret = 0;

  for(int i = 1; i <= numTS; i++)
    {
      x.Rows(1,sizeTS) = p_ts.Column(i);
      FFT(x, dummy, fft_real, fft_im);
      
      for(int j = 1; j <= p_zeropad; j++)
	{
	  // (x+iy)(x-iy) = x^2 + y^2
	  fft_real(j) = fft_real(j)*fft_real(j) + fft_im(j)*fft_im(j);
	  fft_im(j) = 0;
	}
      
      FFTI(fft_real, fft_im, realifft, dummy2);
      
      float varx = var(x.Rows(1,sizeTS)).AsScalar();
      ret.Column(i) = realifft.Rows(1,lag);
      
      for(int j = 1; j <= lag-1; j++)
	{
	  // Correction to make autocorr unbiased and normalised
	  ret(j,i) = ret(j,i)/((sizeTS-j)*varx);
	}  
    }
}

ReturnMatrix xcorr(const Matrix& p_ts, int lag, int p_zeropad )
{
  Matrix r;
  
  xcorr(p_ts,r,lag,p_zeropad);
  r.Release();
  return r;
}

void detrend(Matrix& p_ts, int p_level)
{
  Tracer trace("MISCMATHS::detrend");

  int sizeTS = p_ts.Nrows();     

  // p_ts = b*a + e (OLS regression)
  // e is detrended data
  Matrix a(sizeTS, p_level+1);  

  // Create a
  for(int t = 1; t <= sizeTS; t++)
    {
      for(int l = 0; l <= p_level; l++)
	a(t,l+1) = pow((float)t/sizeTS,l);
    }
          
  // Form residual forming matrix R:
  Matrix R = IdentityMatrix(sizeTS)-a*pinv(a);

  for(int t = 1; t <= sizeTS; t++)
    {
      p_ts.Column(t) = R*p_ts.Column(t);
    }      
}



ReturnMatrix read_vest(string p_fname)
{
  ifstream in;
  in.open(p_fname.c_str(), ios::in);
  
  if(!in) throw Exception(string("Unable to open "+p_fname).c_str());
  
  int numWaves = 0;
  int numPoints = 0;
  
  string str;
  
  while(true)
    {
      if(!in.good()) throw Exception(string(p_fname+" is not a valid vest file").c_str());	
      in >> str;
      if(str == "/Matrix")
	break;
      else if(str == "/NumWaves")
	{
	  in >> numWaves;
	}
      else if(str == "/NumPoints" || str == "/NumContrasts")
	{
	  in >> numPoints;
	}
    }
  
  Matrix p_mat(numPoints, numWaves);
  
  for(int i = 1; i <= numPoints; i++)
    {
      for(int j = 1; j <= numWaves; j++)    
	{
	  if (!in.eof()) in >> ws >> p_mat(i,j) >> ws;
	  else throw Exception(string(p_fname+" has insufficient data points").c_str());
	}
    }
  
  in.close();

  p_mat.Release();
  return p_mat;
}

void ols(const Matrix& data,const Matrix& des,const Matrix& tc, Matrix& cope,Matrix& varcope){
  // ols
  // data is t x v
  // des is t x ev (design matrix)
  // tc is cons x ev (contrast matrix)
  // cope and varcope will be cons x v
  // but will be resized if they are wrong
  // hence may be passed in uninitialised
  // TB 2004
  if(data.Nrows() != des.Nrows()){
    cerr <<"MISCMATHS::ols - data and design have different number of time points"<<endl;
    exit(-1);
  }
  if(des.Ncols() != tc.Ncols()){
    cerr <<"MISCMATHS::ols - design and contrast matrix have different number of EVs"<<endl;
    exit(-1);
  }  
  Matrix pdes = pinv(des);
  Matrix prevar=diag(tc*pdes*pdes.t()*tc.t());
  Matrix R=IdentityMatrix(des.Nrows())-des*pdes;
  float tR=R.Trace();
  Matrix pe=pdes*data;
  cope=tc*pe;
  Matrix res=data-des*pe;
  Matrix sigsq=sum(SP(res,res))/tR;
  varcope=prevar*sigsq;
  
  
}

float ols_dof(const Matrix& des){
  Matrix pdes = pinv(des);
  Matrix R=IdentityMatrix(des.Nrows())-des*pdes;
  return R.Trace();
}

int conjgrad(ColumnVector& x, const Matrix& A, const ColumnVector& b, int maxit,
	     float reltol)
{
  // solves:  A * x = b    (for x)
  // implementation of algorithm in Golub and Van Loan (3rd ed, page 527)
  ColumnVector rk1, rk2, pk, apk;
  double betak, alphak, rk1rk1=0, rk2rk2, r00=0;
  int k=0;
  rk1 = b - A*x;  // a *big* calculation
  for (int n=1; n<=maxit; n++) {
    k++;
    if (k==1) {
      pk = rk1;
      rk1rk1 = (rk1.t() * rk1).AsScalar();
      r00=rk1rk1;
    } else {
      rk2rk2 = rk1rk1;  // from before
      rk1rk1 = (rk1.t() * rk1).AsScalar();
      if (rk2rk2<1e-10*rk1rk1) {
	cerr << "WARNING:: Conj Grad - low demoninator (rk2rk2)" << endl; 
	if (rk2rk2<=0) {
	  cerr << "Aborting conj grad ..." << endl;
	  return 1;  
	}
      }
      betak = rk1rk1 / rk2rk2;
      pk = rk1 + betak * pk;  // note RHS pk is p(k-1) in algorithm
    }  
    // stop if sufficient accuracy is achieved
    if (rk1rk1<reltol*reltol*r00) return 0;

    apk = A * pk;  // the *big* calculation in this algorithm

    ColumnVector pap = pk.t() * apk;
    if (pap.AsScalar()<0) { 
      cerr << "ERROR:: Conj Grad - negative eigenvector found (matrix must be symmetric and positive-definite)\nAborting ... " << endl; 
      return 2;
    } else if (pap.AsScalar()<1e-10) { 
      cerr << "WARNING:: Conj Grad - nearly null eigenvector found (terminating early)" << endl; 
      return 1;
    } else {
      alphak = rk1rk1 / pap.AsScalar();
    }
    x = x + alphak * pk;  // note LHS is x(k) and RHS is x(k-1) in algorithm
    rk2 = rk1;  // update prior to the next step
    rk1 = rk1 - alphak * apk;  // note LHS is r(k) in algorithm
  }
  return 0;
}

int conjgrad(ColumnVector& x, const Matrix& A, const ColumnVector& b, int maxit)
{ 
  return conjgrad(x,A,b,maxit,1e-10);
}

float csevl(const float x, const ColumnVector& cs, const int n)
  {
 
    float b0 = 0;
    float b1 = 0;
    float b2 = 0;
    const float twox=2*x;
    
    for(int i=1; i<=n; i++)
      {
	b2=b1;
	b1=b0;
	b0=twox*b1-b2+cs(n+1-i);
      }
    
    return 0.5*(b0-b2);
  }

  float digamma(const float x)
  { 
    int ntapsi(16);
    int ntpsi(23);
    ColumnVector psics(ntpsi);
    ColumnVector apsics(ntapsi);

    psics << -.038057080835217922E0<<
      .49141539302938713E0<<
      -.056815747821244730E0<<
      .008357821225914313E0<<
      -.001333232857994342E0<<
      .000220313287069308E0<<
      -.000037040238178456E0<<
      .000006283793654854E0<<
      -.000001071263908506E0<<
      .000000183128394654E0<<
      -.000000031353509361E0<<
      .000000005372808776E0<<
      -.000000000921168141E0<<
      .000000000157981265E0<<
      -.000000000027098646E0<<
      .000000000004648722E0<<
      -.000000000000797527E0<<
      .000000000000136827E0<<
      -.000000000000023475E0<<
      .000000000000004027E0<<
      -.000000000000000691E0<<
      .000000000000000118E0<<
      -.000000000000000020E0;
	  
    apsics <<-.0204749044678185E0<<
      -.0101801271534859E0<<
      .0000559718725387E0<<
      -.0000012917176570E0<<
      .0000000572858606E0<<
      -.0000000038213539E0<<
      .0000000003397434E0<<
      -.0000000000374838E0<<
      .0000000000048990E0<<
      -.0000000000007344E0<<
      .0000000000001233E0<<
      -.0000000000000228E0<<
      .0000000000000045E0<<
      -.0000000000000009E0<<
      .0000000000000002E0<<
      -.0000000000000000E0;

    float y = fabs(x);
    float psi;

    if(y<2.0)
      {
	// do we need to deal with the following case?
	// c psi(x) for -2. .lt. x .lt. 2.

	int n = int(floor(x));
	y = x - n;
	n = n - 1;
	psi = csevl(2*y-1, psics, ntpsi);
	if(n==-1)
	  {
	    psi = psi - 1.0/x;
	  }
      }
    else
      {
	const float aux = csevl(8/(Sqr(y))-1, apsics, ntapsi);
	psi = log(x) - 0.5/x + aux;
      }
    
    return psi;
  }

  void glm_vb(const Matrix& X, const ColumnVector& Y, ColumnVector& B, SymmetricMatrix& ilambda_B, int niters)
  {
    // Does Variational Bayes inference on GLM Y=XB+e with ARD priors on B
    // design matrix X should be num_tpts*num_evs

    /////////////////////
    // setup
    OUT("Setup");

    int ntpts=Y.Nrows();
    int nevs=X.Ncols();

    if(ntpts!=X.Nrows())
      throw Exception("COCK");

    OUT(nevs);
    OUT(ntpts);

    ColumnVector gam_m(nevs);
    gam_m=1e10;
    float gam_y;

    ColumnVector lambdaB(nevs);
    if(nevs<ntpts-10)
      {
	// initialise with OLS
	B=pinv(X)*Y;
	ColumnVector res=Y-X*B;
	gam_y=(ntpts-nevs)/(res.t()*res).AsScalar();

	ilambda_B << (X.t()*X*gam_y).i();
	lambdaB=0;
	for(int l=1; l <= nevs; l++)
	  {
	    lambdaB(l)=ilambda_B(l,l);
	  }
      }
    else
      {
	OUT("no ols");
	B.ReSize(nevs);
	B=0;
	lambdaB=1;

// 	ColumnVector res=Y-X*B;
// 	gam_y=ntpts/(res.t()*res).AsScalar();

	gam_y=10;
      }

//     OUT(B(1));
//     OUT(lambdaB(1));

    float trace_ilambdaZZ=1;

    SymmetricMatrix ZZ;
    ZZ << X.t()*X;

    Matrix ZY = X.t()*Y;

    float YY=0;
    for(int t=1; t <= ntpts; t++)
      YY += Sqr(Y(t));

    /////////////////////
    // iterate
    OUT("Iterate");

    int i = 1;;
    for(; i<=niters; i++)
      {
	cout<<i<<",";
	////////////////////
	// update phim
	for(int l=1; l <= nevs; l++)
	  {
	    float b_m0 = 1e10;
	    float c_m0 = 2;

	    float c_m = 1.0/2.0 + c_m0;	    
	    float b_m = 1.0/(0.5*(Sqr(B(l))+lambdaB(l))+1.0/b_m0);
	    gam_m(l) = b_m*c_m;	    
	  }

// 	OUT(gam_m(1));

	////////////////////
	// update B
	ColumnVector beta(nevs);
	beta = 0;
	SymmetricMatrix lambda_B(nevs);
	lambda_B = 0;

	for(int l=1; l <= nevs; l++)
	  lambda_B(l,l)=gam_m(l);

	SymmetricMatrix tmp = lambda_B + gam_y*ZZ;
	lambda_B << tmp;

	beta += gam_y*ZY;

	ilambda_B << lambda_B.i();
	B=ilambda_B*beta;

	lambdaB.ReSize(nevs);
	lambdaB=0;
	for(int l=1; l <= nevs; l++)
	  {
	    lambdaB(l)=ilambda_B(l,l);
	  }
	
	////////////////////
	// compute trace for noise precision phiy update
	
	SymmetricMatrix tmp3;
	tmp3 << ilambda_B;
	
	SymmetricMatrix tmp2;
	tmp2 << tmp3*ZZ;
	
	trace_ilambdaZZ=tmp2.Trace();	
// 	OUT(trace_ilambdaZZ);


	/////////////////////
	// update phiy
	float b_y0 = 1e10;
	float c_y0 = 1;

	float c_y = (ntpts-1)/2.0 + c_y0;

	float sum = YY + (B.t()*ZZ*B).AsScalar() - 2*(B.t()*ZY).AsScalar();
	
	float b_y = 1.0/(0.5*(sum + trace_ilambdaZZ)+1/b_y0);
	
	gam_y = b_y*c_y;

// 	OUT(gam_y);	     

      }

    cout << endl;
  }

vector<float> ColumnVector2vector(const ColumnVector& col)
{
  vector<float> vec(col.Nrows());
  
  for(int c = 0; c < col.Nrows(); c++)
    vec[c] = col(c+1);
  
  return vec;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

// Uninteresting byte swapping functions

typedef struct { unsigned char a,b ; } TWObytes ;

void Swap_2bytes( int n , void *ar )    /* 2 bytes at a time */
{
  register int ii ;
   register TWObytes *tb = (TWObytes *)ar ;
   register unsigned char tt ;

   for( ii=0 ; ii < n ; ii++ ){
     tt = tb[ii].a ; tb[ii].a = tb[ii].b ; tb[ii].b = tt ;
   }
   return ;
}

/*---------------------------------------------------------------------------*/

typedef struct { unsigned char a,b,c,d ; } FOURbytes ;

void Swap_4bytes( int n , void *ar )    /* 4 bytes at a time */
{
   register int ii ;
   register FOURbytes *tb = (FOURbytes *)ar ;
   register unsigned char tt ;

   for( ii=0 ; ii < n ; ii++ ){
     tt = tb[ii].a ; tb[ii].a = tb[ii].d ; tb[ii].d = tt ;
     tt = tb[ii].b ; tb[ii].b = tb[ii].c ; tb[ii].c = tt ;
   }
   return ;
}

/*---------------------------------------------------------------------------*/

typedef struct { unsigned char a,b,c,d , D,C,B,A ; } EIGHTbytes ;

void Swap_8bytes( int n , void *ar )    /* 8 bytes at a time */
{
   register int ii ;
   register EIGHTbytes *tb = (EIGHTbytes *)ar ;
   register unsigned char tt ;

   for( ii=0 ; ii < n ; ii++ ){
     tt = tb[ii].a ; tb[ii].a = tb[ii].A ; tb[ii].A = tt ;
     tt = tb[ii].b ; tb[ii].b = tb[ii].B ; tb[ii].B = tt ;
     tt = tb[ii].c ; tb[ii].c = tb[ii].C ; tb[ii].C = tt ;
     tt = tb[ii].d ; tb[ii].d = tb[ii].D ; tb[ii].D = tt ;
   }
   return ;
}

/*---------------------------------------------------------------------------*/

typedef struct { unsigned char a,b,c,d,e,f,g,h ,
                               H,G,F,E,D,C,B,A  ; } SIXTEENbytes ;

void Swap_16bytes( int n , void *ar )    /* 16 bytes at a time */
{
   register int ii ;
   register SIXTEENbytes *tb = (SIXTEENbytes *)ar ;
   register unsigned char tt ;

   for( ii=0 ; ii < n ; ii++ ){
     tt = tb[ii].a ; tb[ii].a = tb[ii].A ; tb[ii].A = tt ;
     tt = tb[ii].b ; tb[ii].b = tb[ii].B ; tb[ii].B = tt ;
     tt = tb[ii].c ; tb[ii].c = tb[ii].C ; tb[ii].C = tt ;
     tt = tb[ii].d ; tb[ii].d = tb[ii].D ; tb[ii].D = tt ;

     tt = tb[ii].e ; tb[ii].e = tb[ii].E ; tb[ii].E = tt ;
     tt = tb[ii].f ; tb[ii].f = tb[ii].F ; tb[ii].F = tt ;
     tt = tb[ii].g ; tb[ii].g = tb[ii].G ; tb[ii].G = tt ;
     tt = tb[ii].h ; tb[ii].h = tb[ii].H ; tb[ii].H = tt ;
   }
   return ;
}

/*---------------------------------------------------------------------------*/

void Swap_Nbytes( int n , int siz , void *ar )  /* subsuming case */
{
   switch( siz ){
     case 2:  Swap_2bytes ( n , ar ) ; break ;
     case 4:  Swap_4bytes ( n , ar ) ; break ;
     case 8:  Swap_8bytes ( n , ar ) ; break ;
     case 16: Swap_16bytes( n , ar ) ; break ;
   }
   return ;
}

// end namespace MISCMATHS
}
