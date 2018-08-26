/*  General IO functions (images and transformation files)

    Mark Jenkinson and Matthew Webster, FMRIB Image Analysis Group

    Copyright (C) 2000-2008 University of Oxford  */

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


#if !defined(__newimageio_h)
#define __newimageio_h

#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include "newmatio.h"
#include "newimage.h"
#include "complexvolume.h"
#include "fslio/fslio.h"
#include "miscmaths/miscmaths.h"

using namespace NEWMAT;

namespace NEWIMAGE {

string fslbasename(const string& filename);
int make_basename(string& filename);
int find_pathname(string& filename);
bool fsl_imageexists(const string& filename);

  // read

template <class T>
int read_volume(volume<T>& target, const string& filename);
template <class T>
int read_volumeROI(volume<T>& target, const string& filename,
		   int x0, int y0, int z0, int x1, int y1, int z1);
template <class T>
int read_volumeROI(volume<T>& target, const string& filename,
		   int x0, int y0, int z0, int x1, int y1, int z1, 
		   int xskip, int yskip, int zskip);
template <class T>
int read_volume_hdr_only(volume<T>& target, const string& filename);

int read_complexvolume(volume<float>& realvol, volume<float>& imagvol,
		       const string& filename);
int read_complexvolume(complexvolume& vol, const string& filename);


template <class T>
int read_volume4D(volume4D<T>& target, const string& filename);
template <class T>
int read_volume4DROI(volume4D<T>& target, const string& filename, 
		     int x0, int y0, int z0, int t0, 
		     int x1, int y1, int z1, int t1);
template <class T>
int read_volume4DROI(volume4D<T>& target, const string& filename, 
		     int x0, int y0, int z0, int t0, 
		     int x1, int y1, int z1, int t1,
		     int xskip, int yskip, int zskip, int tskip);
template <class T>
int read_volume4D_hdr_only(volume4D<T>& target, const string& filename);
int read_complexvolume4D(volume4D<float>& realvol, volume4D<float>& imagvol,
			 const string& filename);
int read_complexvolume4D(complexvolume& vol, const string& filename);


  // save

template <class T>
int save_volume(const volume<T>& source, const string& filename);
int save_complexvolume(const volume<float>& realvol, 
		       const volume<float>& imagvol, const string& filename);
int save_complexvolume(const complexvolume& vol, const string& filename);


template <class T>
int save_volume4D(const volume4D<T>& source, const string& filename);
int save_complexvolume4D(const volume4D<float>& realvol, 
			 const volume4D<float>& imagvol, 
			 const string& filename);
int save_complexvolume4D(const complexvolume& vol, const string& filename);

template <class T>
int save_volume_datatype(const volume<T>& source, const string& filename,
			 short datatype);
template <class T>
int save_volume4D_datatype(const volume4D<T>& source, const string& filename,
			   short datatype);

template <class T>
int save_volume_filetype(const volume<T>& source, const string& filename,
			 int filetype);
template <class T>
int save_volume4D_filetype(const volume4D<T>& source, const string& filename,
			   int filetype);

// Helper functions
short closestTemplatedType(const short inputType);

short dtype(const char* T);
short dtype(const short* T);
short dtype(const int* T);
short dtype(const float* T);
short dtype(const double* T);

short dtype(const volume<char>& vol);
short dtype(const volume<short>& vol);
short dtype(const volume<int>& vol);
short dtype(const volume<float>& vol);
short dtype(const volume<double>& vol);

short dtype(const volume4D<char>& vol);
short dtype(const volume4D<short>& vol);
short dtype(const volume4D<int>& vol);
short dtype(const volume4D<float>& vol);
short dtype(const volume4D<double>& vol);

short dtype(const string& filename);

short fslFileType(const string& filename);

// Boring overloads to enable different names (load and write)


// load

template <class T>
int load_volume(volume<T>& target, const string& filename);
template <class T>
int load_volume_hdr_only(volume<T>& target, const string& filename);
template <class T>
int load_complexvolume(volume<float>& realvol, volume<float>& imagvol,
		       const string& filename);
int load_complexvolume(complexvolume& vol, const string& filename);

template <class T>
int load_volume4D(volume4D<T>& target, const string& filename);
template <class T>
int load_volume4D_hdr_only(volume4D<T>& target, const string& filename);
int load_complexvolume4D(volume4D<float>& realvol, volume4D<float>& imagvol,
			 const string& filename);
// write

template <class T>
int write_volume(const volume<T>& source, const string& filename);
int write_complexvolume(const volume<float>& realvol, 
		       const volume<float>& imagvol, const string& filename);
int write_complexvolume(const complexvolume& vol, const string& filename);

template <class T>
int write_volume4D(const volume4D<T>& source, const string& filename);
int write_complexvolume4D(const volume4D<float>& realvol, 
			 const volume4D<float>& imagvol, 
			 const string& filename);

template <class T>
int write_volume_datatype(const volume<T>& source, const string& filename,
			 short datatype);
template <class T>
int write_volume4D_datatype(const volume4D<T>& source, const string& filename,
			   short datatype);
template <class T>
int write_volume_filetype(const volume<T>& source, const string& filename,
			 int filetype);
template <class T>
int write_volume4D_filetype(const volume4D<T>& source, const string& filename,
			   int filetype);




////////////////////////////////////////////////////////////////////////
///////////////////////// TEMPLATE DEFINITIONS /////////////////////////
////////////////////////////////////////////////////////////////////////

template <class T>
void FslReadBuffer(FSLIO* IP, T* tbuffer) 
{
  short sx,sy,sz,st;
  FslGetDim(IP,&sx,&sy,&sz,&st);
  size_t imagesize=sx*sy*sz;
  short type;
  float slope, intercept;
  bool doscaling = false;

  FslGetDataType(IP,&type);
  doscaling = FslGetIntensityScaling(IP,&slope,&intercept);
  if ( (dtype(tbuffer) == type) && (!doscaling) ) {
    FslReadVolumes(IP,tbuffer,1);
  } else {
    switch(type)
      {
      case DT_SIGNED_SHORT:
	{
	  short* sbuffer;
	  try {
	    sbuffer=new short[imagesize];
	  } catch(...) { sbuffer=0; }
	  if (sbuffer==0) { imthrow("Out of memory for size "+num2str(imagesize),99); }
	  FslReadVolumes(IP,sbuffer,1);
	  if (doscaling) convertbuffer(sbuffer,tbuffer,imagesize,slope,intercept);
	  else convertbuffer(sbuffer,tbuffer,imagesize);
	  delete[] sbuffer;
	}
	break;
      case DT_UNSIGNED_CHAR:
	{
	  unsigned char* sbuffer;
	  try {
	    sbuffer=new unsigned char[imagesize];
	  } catch(...) { sbuffer=0; }
	  if (sbuffer==0) { imthrow("Out of memory for size "+num2str(imagesize),99); }
	  FslReadVolumes(IP,sbuffer,1);
	  if (doscaling) convertbuffer(sbuffer,tbuffer,imagesize,slope,intercept);
	  else convertbuffer(sbuffer,tbuffer,imagesize);
	  delete[] sbuffer;
	}
	break;
      case DT_SIGNED_INT:
	{
	  int* sbuffer;
	  try {
	    sbuffer=new int[imagesize];
	  } catch(...) { sbuffer=0; }
	  if (sbuffer==0) { imthrow("Out of memory for size "+num2str(imagesize),99); }
	  FslReadVolumes(IP,sbuffer,1);
	  if (doscaling) convertbuffer(sbuffer,tbuffer,imagesize,slope,intercept);
	  else convertbuffer(sbuffer,tbuffer,imagesize);
	  delete[] sbuffer;
	}
	break;
      case DT_FLOAT:
	{
	  float* sbuffer;
	  try {
	    sbuffer=new float[imagesize];
	  } catch(...) { sbuffer=0; }
	  if (sbuffer==0) { imthrow("Out of memory for size "+num2str(imagesize),99); }
	  FslReadVolumes(IP,sbuffer,1);
	  if (doscaling) convertbuffer(sbuffer,tbuffer,imagesize,slope,intercept);
	  else convertbuffer(sbuffer,tbuffer,imagesize);
	  delete[] sbuffer;
	}
	break;
      case DT_DOUBLE:
	{
	  double* sbuffer;
	  try {
	    sbuffer=new double[imagesize];
	  } catch(...) { sbuffer=0; }
	  if (sbuffer==0) { imthrow("Out of memory for size "+num2str(imagesize),99); }
	  FslReadVolumes(IP,sbuffer,1);
	  if (doscaling) convertbuffer(sbuffer,tbuffer,imagesize,slope,intercept);
	  else convertbuffer(sbuffer,tbuffer,imagesize);
	  delete[] sbuffer;
	}
	break;
	/*------------------- new codes for NIFTI ---*/
      case DT_INT8:
	{
	  signed char* sbuffer;
	  try {
	    sbuffer=new signed char[imagesize];
	  } catch(...) { sbuffer=0; }
	  if (sbuffer==0) { imthrow("Out of memory for size "+num2str(imagesize),99); }
	  FslReadVolumes(IP,sbuffer,1);
	  if (doscaling) convertbuffer(sbuffer,tbuffer,imagesize,slope,intercept);
	  else convertbuffer(sbuffer,tbuffer,imagesize);
	  delete[] sbuffer;
	}
	break;
      case DT_UINT16:
	{
	  unsigned short* sbuffer;
	  try {
	    sbuffer=new unsigned short[imagesize];
	  } catch(...) { sbuffer=0; }
	  if (sbuffer==0) { imthrow("Out of memory for size "+num2str(imagesize),99); }
	  FslReadVolumes(IP,sbuffer,1);
	  if (doscaling) convertbuffer(sbuffer,tbuffer,imagesize,slope,intercept);
	  else convertbuffer(sbuffer,tbuffer,imagesize);
	  delete[] sbuffer;
	}
	break;
      case DT_UINT32:
	{
	  unsigned int* sbuffer;
	  try {
	    sbuffer=new unsigned int[imagesize];
	  } catch(...) { sbuffer=0; }
	  if (sbuffer==0) { imthrow("Out of memory for size "+num2str(imagesize),99); }
	  FslReadVolumes(IP,sbuffer,1);
	  if (doscaling) convertbuffer(sbuffer,tbuffer,imagesize,slope,intercept);
	  else convertbuffer(sbuffer,tbuffer,imagesize);
	  delete[] sbuffer;
	}
	break;
      case DT_INT64:
	{
	  long signed int* sbuffer;
	  try {
	    sbuffer=new long signed int[imagesize];
	  } catch(...) { sbuffer=0; }
	  if (sbuffer==0) { imthrow("Out of memory for size "+num2str(imagesize),99); }
	  FslReadVolumes(IP,sbuffer,1);
	  if (doscaling) convertbuffer(sbuffer,tbuffer,imagesize,slope,intercept);
	  else convertbuffer(sbuffer,tbuffer,imagesize);
	  delete[] sbuffer;
	}
	break;
      case DT_UINT64:
	{
	  long unsigned int* sbuffer;
	  try {
	    sbuffer=new long unsigned int[imagesize];
	  } catch(...) { sbuffer=0; }
	  if (sbuffer==0) { imthrow("Out of memory for size "+num2str(imagesize),99); }
	  FslReadVolumes(IP,sbuffer,1);
	  if (doscaling) convertbuffer(sbuffer,tbuffer,imagesize,slope,intercept);
	  else convertbuffer(sbuffer,tbuffer,imagesize);
	  delete[] sbuffer;
	}
	break;
      default:
	  /* includes: DT_BINARY, DT_RGB, DT_ALL, DT_FLOAT128, DT_COMPLEX's */
	ostringstream errmsg;
	errmsg << "Fslread: DT " << type <<  " not supported";
	perror(errmsg.str().c_str());
      }
  }
}


void check_filename(const string& basename);

FSLIO* NewFslOpen(const string& filename, const string& permissions, 
		  int filetype);

FSLIO* NewFslOpen(const string& filename, const string& permissions);

// External functions

// READ FUNCTIONS


template <class T>
  void set_volume_properties(FSLIO* IP1, volume<T>& target);



template <class T>
int read_volumeROI(volume<T>& target, const string& filename, 
		   int x0, int y0, int z0, int x1, int y1, int z1,
		   int xskip, int yskip, int zskip)
{
  int retval=read_volumeROI(target,filename,x0,y0,z0,x1,y1,z1);
  if (retval==0) {
    if (xskip<1) xskip=1;    
    if (yskip<1) yskip=1;
    if (zskip<1) zskip=1;
    int sx=(target.maxx()-target.minx())/xskip + 1;
    int sy=(target.maxy()-target.miny())/yskip + 1;
    int sz=(target.maxz()-target.minz())/zskip + 1;
    volume<T> tmpvol(sx,sy,sz);
    int xx=0, yy=0, zz=0, x=0, y=0, z=0;
    for (z=target.minz(), zz=0; z<=target.maxz(); z+=zskip, zz++) {
      for (y=target.miny(), yy=0; y<=target.maxy(); y+=yskip, yy++) {
	for (x=target.minx(), xx=0; x<=target.maxx(); x+=xskip, xx++) {
	  tmpvol(xx,yy,zz) = target(x,y,z);
	}
      }
    }
    tmpvol.copyproperties(target);
    target = tmpvol;
  }
  return retval;
}




template <class T>
int read_volumeROI(volume<T>& target, const string& filename, 
		   short& dtype, bool read_img_data,
		   int x0, int y0, int z0, int x1, int y1, int z1,
		   bool swap2radiological=true);


template <class T>
int read_volumeROI(volume<T>& target, const string& filename, 
		     int x0, int y0, int z0, 
		     int x1, int y1, int z1)
{
  short dtype;
  return read_volumeROI(target,filename,dtype,true,
			  x0,y0,z0,x1,y1,z1);
}


template <class T>
int read_volume(volume<T>& target, const string& filename,short& dtype, bool read_img_data)
{
  return read_volumeROI(target,filename,dtype,read_img_data,0,0,0,-1,-1,-1);
}

template <class T>
int read_volume(volume<T>& target, const string& filename)
{
  short dtype;
  int retval = read_volume(target,filename,dtype,true);
  return retval;
}

template <class T>
int read_volume_hdr_only(volume<T>& target, const string& filename)
{
  short dtype;
  int retval = read_volume(target,filename,dtype,false);
  return retval;
}

template <class T>
int read_volume4DROI(volume4D<T>& target, const string& filename, 
		     int x0, int y0, int z0, int t0, 
		     int x1, int y1, int z1, int t1,
		     int xskip, int yskip, int zskip, int tskip)
{
  int retval=read_volume4DROI(target,filename,x0,y0,z0,t0,x1,y1,z1,t1);
  if (retval==0) {
    if (xskip<1) xskip=1;    
    if (yskip<1) yskip=1;
    if (zskip<1) zskip=1;
    if (tskip<1) tskip=1;
    int sx=(target.maxx()-target.minx())/xskip + 1;
    int sy=(target.maxy()-target.miny())/yskip + 1;
    int sz=(target.maxz()-target.minz())/zskip + 1;
    int st=(target.maxt()-target.mint())/tskip + 1;
    volume4D<T> tmpvol(sx,sy,sz,st);
    int xx=0, yy=0, zz=0, tt=0, x=0, y=0, z=0, t=0;
    for (t=target.mint(), tt=0; t<=target.maxt(); t+=tskip, tt++) {
      for (z=target.minz(), zz=0; z<=target.maxz(); z+=zskip, zz++) {
	for (y=target.miny(), yy=0; y<=target.maxy(); y+=yskip, yy++) {
	  for (x=target.minx(), xx=0; x<=target.maxx(); x+=xskip, xx++) {
	    tmpvol(xx,yy,zz,tt) = target(x,y,z,t);
	  }
	}
      }
    }
    tmpvol.copyproperties(target[0]);
    target = tmpvol;
  }
  return retval;
}


template <class T>
int read_volume4DROI(volume4D<T>& target, const string& filename, 
		     short& dtype, bool read_img_data,
		     int x0, int y0, int z0, int t0, 
		     int x1, int y1, int z1, int t1,
		     bool swap2radiological=true);

template <class T>
int read_volume4DROI(volume4D<T>& target, const string& filename, 
		     int x0, int y0, int z0, int t0, 
		     int x1, int y1, int z1, int t1)
{
  short dtype;
  return read_volume4DROI(target,filename,dtype,true,
			  x0,y0,z0,t0,x1,y1,z1,t1);
}

template <class T>
int read_volume4D(volume4D<T>& target, const string& filename, 
		  short& dtype, bool read_img_data)
{
  return read_volume4DROI(target,filename,dtype,read_img_data,
			  0,0,0,0,-1,-1,-1,-1);
}


template <class T>
int read_volume4D(volume4D<T>& target, const string& filename)
{
  short dtype;
  int retval = read_volume4D(target,filename,dtype,true);
  return retval;
}


template <class T>
int read_volume4D_hdr_only(volume4D<T>& target, const string& filename)
{
  short dtype;
  int retval = read_volume4D(target,filename,dtype,false);
  return retval;
}

// SAVE FUNCTIONS


mat44 newmat2mat44(const Matrix& nmat);


template <class T>
int set_fsl_hdr(const volume<T>& source, FSLIO *OP, int tsize, float tdim, float scalingSlope=1.0) //temp bool while converting
{
  Tracer tr("set_fsl_hdr");
    
  FslSetDim(OP,source.xsize(),source.ysize(),source.zsize(),tsize);
  FslSetDataType(OP, dtype(source));
  FslSetVoxDim(OP,source.xdim(), source.ydim(), source.zdim(), tdim);

  FslSetStdXform(OP,source.sform_code(),newmat2mat44(source.sform_mat()));
  FslSetRigidXform(OP,source.qform_code(),newmat2mat44(source.qform_mat()));
  
  FslSetIntent(OP,source.intent_code(),source.intent_param(1),
	       source.intent_param(2),source.intent_param(3));
  FslSetIntensityScaling(OP,scalingSlope,0.0);
  FslSetCalMinMax(OP,source.getDisplayMinimum(),source.getDisplayMaximum());
  FslSetAuxFile(OP,source.getAuxFile().c_str());
  return 0;
}

template <class T>
int save_basic_volume(const volume<T>& source, const string& filename, 
		      int filetype, bool save_orig=false);

template <class T>
int save_basic_volume4D(const volume4D<T>& source, const string& filename,
			int filetype, bool save_orig=false);

template <class T>
int save_volume(const volume<T>& source, const string& filename)
{
  return save_basic_volume(source,fslbasename(filename),-1);
}
template <class T>
int save_volume4D(const volume4D<T>& source, const string& filename)
{
  return save_basic_volume4D(source,fslbasename(filename),-1);
}


template <class T>
int save_volume_datatype(const volume<T>& source, const string& filename,short datatype)
{
  datatype=closestTemplatedType(datatype);
  if (dtype(source) == datatype) {
    return save_volume(source,filename);
  } else {
    switch(datatype)
      {
      case DT_SIGNED_SHORT:
	{
	  volume<short> svol;
	  copyconvert(source,svol);
	  return save_volume(svol,filename);
	}
	break;
      case DT_UNSIGNED_CHAR:
	{
	  volume<char> svol;
	  copyconvert(source,svol);
	  return save_volume(svol,filename);
	}
	break;
      case DT_SIGNED_INT:
	{
	  volume<int> svol;
	  copyconvert(source,svol);
	  return save_volume(svol,filename);
	}
	break;
      case DT_FLOAT:
	{
	  volume<float> svol;
	  copyconvert(source,svol);
	  return save_volume(svol,filename);
	}
	break;
      case DT_DOUBLE:
	{
	  volume<double> svol;
	  copyconvert(source,svol);
	  return save_volume(svol,filename);
	}
	break;
      default:
	ostringstream errmsg;
	errmsg << "Fslread: DT " << datatype <<  " not supported";
	perror(errmsg.str().c_str());
      }
  }
  return -1;  // should never get here
}
  

template <class T>
int save_volume4D_datatype(const volume4D<T>& source, const string& filename,short datatype)
{
  datatype=closestTemplatedType(datatype);
  if (dtype(source) == datatype) {
    return save_volume4D(source,filename);
  } else {
    switch(datatype)
      {
      case DT_SIGNED_SHORT:
	{
	  volume4D<short> svol;
	  copyconvert(source,svol);
	  return save_volume4D(svol,filename);
	}
	break;
      case DT_UNSIGNED_CHAR:
	{
	  volume4D<char> svol;
	  copyconvert(source,svol);
	  return save_volume4D(svol,filename);
	}
	break;
      case DT_SIGNED_INT:
	{
	  volume4D<int> svol;
	  copyconvert(source,svol);
	  return save_volume4D(svol,filename);
	}
	break;
      case DT_FLOAT:
	{
	  volume4D<float> svol;
	  copyconvert(source,svol);
	  return save_volume4D(svol,filename);
	}
	break;
      case DT_DOUBLE:
	{
	  volume4D<double> svol;
	  copyconvert(source,svol);
	  return save_volume4D(svol,filename);
	}
	break;
      default:
	ostringstream errmsg;
	errmsg << "Fslread: DT " << datatype <<  " not supported";
	perror(errmsg.str().c_str());
      }
  }
  return -1;  // should never get here
}

// old versions call _dtype - kept for compatability
template <class T>
int save_volume_dtype(const volume<T>& source, const string& filename,
			 short datatype)
{
  return save_volume_datatype(source,filename,datatype);
}

template <class T>
int save_volume4D_dtype(const volume4D<T>& source, const string& filename,
			   short datatype)
{
  return save_volume4D_datatype(source,filename,datatype);
}

template <class T>
int save_volume_filetype(const volume<T>& source, const string& filename,
			 int filetype)
{
  return save_basic_volume(source,filename,filetype);
}

template <class T>
int save_volume4D_filetype(const volume4D<T>& source, const string& filename,
			   int filetype)
{
  return save_basic_volume4D(source,filename,filetype);
}

// functions to save without doing any swapping (i.e. just as passed in)

template <class T>
int save_orig_volume(const volume<T>& source, const string& filename)
{
  return save_basic_volume(source,filename,-1,true);
}

template <class T>
int save_orig_volume4D(const volume4D<T>& source, const string& filename)
{
  return save_basic_volume4D(source,filename,-1,true);
}


////////////////////////////////////////////////////////////////////////
///// Boring overloads to enable different names (load and write) //////
////////////////////////////////////////////////////////////////////////

// load

template <class T>
int load_volume(volume<T>& target, const string& filename)
  { return read_volume(target,filename); }

template <class T>
int load_volume_hdr_only(volume<T>& target, const string& filename)
  { return read_volume_hdr_only(target,filename); }

template <class T>
int load_volume4D(volume4D<T>& target, const string& filename)
  { return read_volume4D(target,filename); }

template <class T>
int load_volume4D_hdr_only(volume4D<T>& target, const string& filename)
  { return read_volume4D_hdr_only(target,filename); }

  // write

template <class T>
int write_volume(const volume<T>& source, const string& filename)
  { return save_volume(source,filename); }
template <class T>
int write_volume4D(const volume4D<T>& source, const string& filename)
  { return save_volume4D(source,filename); }
template <class T>
int write_volume_datatype(const volume<T>& source, const string& filename,
			 short datatype)
  { return save_volume_datatype(source,filename,datatype); }
template <class T>
int write_volume4D_datatype(const volume4D<T>& source, const string& filename,
			   short datatype)
  { return save_volume4D_datatype(source,filename,datatype); }
template <class T>
int write_volume_filetype(const volume<T>& source, const string& filename,
			 int filetype)
  { return save_volume_filetype(source,filename,filetype); }
template <class T>
int write_volume4D_filetype(const volume4D<T>& source, const string& filename,
			   int filetype)
  { return save_volume4D_filetype(source,filename,filetype); }


  // Basic I/O functions 
    // read original storage order - do not swap to radiological

template <class T>
int read_orig_volume(volume<T>& target, const string& filename);
template <class T>
int read_orig_volume4D(volume4D<T>& target, const string& filename);
template <class T>
int load_orig_volume(volume<T>& target, const string& filename)
{ return read_orig_volume(target,filename); }
 	 
template <class T>
int load_orig_volume(volume<T>& target, const string& filename);
template <class T>
int load_orig_volume4D(volume4D<T>& target, const string& filename)
{ return read_orig_volume4D(target,filename); }
 	 
 	 
template <class T>
int read_orig_volume(volume<T>& target, const string& filename)
{
  short dtype;
  return read_volumeROI(target,filename,dtype,true,
			0,0,0,-1,-1,-1,false);
}
template <class T>
int read_orig_volume4D(volume4D<T>& target, const string& filename)
{
  short dtype;
  return read_volume4DROI(target,filename,dtype,true,0,0,0,0,-1,-1,-1,-1,false);
}

 	 
}

#endif

