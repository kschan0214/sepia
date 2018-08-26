/*  cluster.cc

    Mark Jenkinson & Matthew Webster, FMRIB Image Analysis Group

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

#ifndef EXPOSE_TREACHEROUS
#define EXPOSE_TREACHEROUS
#endif    

#include <vector>
#include <algorithm>
#include <iomanip>
// #include "newimage/fmribmain.h"
#include "newimage/newimageall.h"
#include "utils/options.h"
// #include "infer.h"
// #include "warpfns/warpfns.h"
// #include "warpfns/fnirt_file_reader.h"

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1        

using namespace NEWIMAGE;
using std::vector;
using namespace Utilities;
using namespace NEWMAT;
using namespace MISCMATHS;


int num(const char x) { return (int) x; }
short int num(const short int x) { return x; }
int num(const int x) { return x; }
float num(const float x) { return x; }
double num(const double x) { return x; }


template <class T>
struct triple { T x; T y; T z; };

template <class T, class S>
void copyconvert(const vector<triple<T> >& oldcoords, 
     vector<triple<S> >& newcoords)
{
  newcoords.erase(newcoords.begin(),newcoords.end());
  newcoords.resize(oldcoords.size());
  for (unsigned int n=0; n<oldcoords.size(); n++) {
    newcoords[n].x = (S) oldcoords[n].x;
    newcoords[n].y = (S) oldcoords[n].y;
    newcoords[n].z = (S) oldcoords[n].z;
  }
}

template <class T>
void MultiplyCoordinateVector(vector<triple<T> >& coords, const Matrix& mat)
{
  ColumnVector vec(4);
  for (unsigned int n=0; n<coords.size(); n++) {
    vec << coords[n].x << coords[n].y << coords[n].z << 1.0;
    vec = mat * vec;     // apply voxel xfm
    coords[n].x = vec(1);
    coords[n].y = vec(2);
    coords[n].z = vec(3);
  }
}

template <class T, class S>
void TransformToReference(vector<triple<T> >& coordlist, const Matrix& affine, 
        const volume<S>& source, const volume<S>& dest, const volume4D<float>& warp,bool doAffineTransform, bool doWarpfieldTransform)
{
  ColumnVector coord(4);
  for (unsigned int n=0; n<coordlist.size(); n++) {
    coord << coordlist[n].x << coordlist[n].y << coordlist[n].z << 1.0;
    if ( doAffineTransform && doWarpfieldTransform ) coord = NewimageCoord2NewimageCoord(affine,warp,true,source,dest,coord);
    if ( doAffineTransform && !doWarpfieldTransform) coord = NewimageCoord2NewimageCoord(affine,source,dest,coord);
    if ( !doAffineTransform && doWarpfieldTransform) coord = NewimageCoord2NewimageCoord(warp,true,source,dest,coord);
    coordlist[n].x = coord(1);
    coordlist[n].y = coord(2);
    coordlist[n].z = coord(3);
  }
}

template <class T>
bool checkIfLocalMaxima(const int& index, const volume<int>& labelim, const volume<T>& zvol, const int& connectivity, const int& x, const int& y, const int& z )
{        
  if (connectivity==6)
    return ( index==labelim(x,y,z) &&
       zvol(x,y,z)>zvol(x,  y,  z-1) &&
       zvol(x,y,z)>zvol(x,  y-1,z) &&
       zvol(x,y,z)>zvol(x-1,y,  z) &&
       zvol(x,y,z)>=zvol(x+1,y,  z) &&
       zvol(x,y,z)>=zvol(x,  y+1,z) &&
       zvol(x,y,z)>=zvol(x,  y,  z+1) );

  else 
    return ( index==labelim(x,y,z) &&
       zvol(x,y,z)>zvol(x-1,y-1,z-1) &&
       zvol(x,y,z)>zvol(x,  y-1,z-1) &&
       zvol(x,y,z)>zvol(x+1,y-1,z-1) &&
       zvol(x,y,z)>zvol(x-1,y,  z-1) &&
       zvol(x,y,z)>zvol(x,  y,  z-1) &&
       zvol(x,y,z)>zvol(x+1,y,  z-1) &&
       zvol(x,y,z)>zvol(x-1,y+1,z-1) &&
       zvol(x,y,z)>zvol(x,  y+1,z-1) &&
       zvol(x,y,z)>zvol(x+1,y+1,z-1) &&
       zvol(x,y,z)>zvol(x-1,y-1,z) &&
       zvol(x,y,z)>zvol(x,  y-1,z) &&
       zvol(x,y,z)>zvol(x+1,y-1,z) &&
       zvol(x,y,z)>zvol(x-1,y,  z) &&
       zvol(x,y,z)>=zvol(x+1,y,  z) &&
       zvol(x,y,z)>=zvol(x-1,y+1,z) &&
       zvol(x,y,z)>=zvol(x,  y+1,z) &&
       zvol(x,y,z)>=zvol(x+1,y+1,z) &&
       zvol(x,y,z)>=zvol(x-1,y-1,z+1) &&
       zvol(x,y,z)>=zvol(x,  y-1,z+1) &&
       zvol(x,y,z)>=zvol(x+1,y-1,z+1) &&
       zvol(x,y,z)>=zvol(x-1,y,  z+1) &&
       zvol(x,y,z)>=zvol(x,  y,  z+1) &&
       zvol(x,y,z)>=zvol(x+1,y,  z+1) &&
       zvol(x,y,z)>=zvol(x-1,y+1,z+1) &&
       zvol(x,y,z)>=zvol(x,  y+1,z+1) &&
       zvol(x,y,z)>=zvol(x+1,y+1,z+1) );

}

template <class T>
void get_stats(const volume<int>& labelim, const volume<T>& origim,
         vector<int>& size,
         vector<T>& maxvals, vector<float>& meanvals,
         vector<triple<int> >& max, vector<triple<float> >& cog,
         bool minv) 
{
  int labelnum = labelim.max();
  size.resize(labelnum+1,0);
  maxvals.resize(labelnum+1, (T) 0);
  meanvals.resize(labelnum+1,0.0f);
  triple<int> zero;
  zero.x = 0; zero.y = 0; zero.z = 0;
  triple<float> zerof;
  zerof.x = 0; zerof.y = 0; zerof.z = 0;
  max.resize(labelnum+1,zero);
  cog.resize(labelnum+1,zerof);
  vector<float> sum(labelnum+1,0.0);
  for (int z=labelim.minz(); z<=labelim.maxz(); z++) {
    for (int y=labelim.miny(); y<=labelim.maxy(); y++) {
      for (int x=labelim.minx(); x<=labelim.maxx(); x++) {
  int idx = labelim(x,y,z);
  T oxyz = origim(x,y,z);
  size[idx]++;
  cog[idx].x+=((float) oxyz)*x;
  cog[idx].y+=((float) oxyz)*y;
  cog[idx].z+=((float) oxyz)*z;
  sum[idx]+=(float) oxyz;
  if ( (size[idx]==1) || 
       ( (oxyz>maxvals[idx]) && (!minv) ) || 
       ( (oxyz<maxvals[idx]) && (minv) ) ) 
    {
      maxvals[idx] = oxyz;
      max[idx].x = x;
      max[idx].y = y;
      max[idx].z = z;
    }
      }
    }
  }
  for (int n=0; n<=labelnum; n++) {
    if (size[n]>0.0) {
      meanvals[n] = (sum[n]/((float) size[n]));
    }
    if (sum[n]>0.0) {
      cog[n].x /= sum[n];
      cog[n].y /= sum[n];
      cog[n].z /= sum[n];
    }
  }
}


vector<int> get_sortindex(const vector<int>& vals)
{
  // return the mapping of old indices to new indices in the
  //   new *ascending* sort of vals
  int length=vals.size();
  vector<pair<int, int> > sortlist(length);
  for (int n=0; n<length; n++) {
    sortlist[n] = pair<int, int>(vals[n],n);
  }
  sort(sortlist.begin(),sortlist.end());  // O(N.log(N))
  vector<int> idx(length);
  for (int n=0; n<length; n++) {
    idx[n] = sortlist[n].second;
  }
  return idx;
}


void get_sizeorder(const vector<int>& size, vector<int>& sizeorder) 
{
  vector<int> sizecopy(size), idx;
  idx = get_sortindex(sizecopy);

  // second part of pair is now the prior-index of the sorted values
  int length = size.size();
  sizeorder.resize(length,0);
  for (int n=0; n<length; n++) {
    sizeorder[idx[n]] = n;    // maps old index to new
  }
}


template <class T, class S>
void relabel_image(const volume<int>& labelim, volume<T>& relabelim,
       const vector<S>& newlabels)
{
  copyconvert(labelim,relabelim);
  for (int z=relabelim.minz(); z<=relabelim.maxz(); z++) 
    for (int y=relabelim.miny(); y<=relabelim.maxy(); y++) 
      for (int x=relabelim.minx(); x<=relabelim.maxx(); x++) 
  relabelim(x,y,z) = (T) newlabels[labelim(x,y,z)];
}


int cluster(int *inputimage, int* outputimage, int* output_nlabel, const int* imgdim, const int num_connect)
{

  volume<int> labelim;
  vector<int> size;
  vector<triple<int> > maxpos;
  vector<triple<float> > cog;
  vector<int> maxvals;
  vector<float> meanvals;

  // read in the volume
  volume<int> mask(imgdim[0], imgdim[1], imgdim[2], inputimage, 0);
  
  // Find the connected components
  labelim = connected_components(mask, num_connect);
  // print_volume_info(labelim,"Labelim");
  
  // process according to the output format
  get_stats(labelim,mask,size,maxvals,meanvals,maxpos,cog,0);
  // cout<<"Number of labels = "<<size.size()<<endl;

  // re-threshold for p
  int length = size.size();
  size[0] = 0;  // force background to be ordered last

  // get sorted index (will revert to cluster size)
  vector<int> idx;
  idx = get_sortindex(size);
  // cout<<size.size()<<" labels in sortedidx"<<endl;

  vector<int> threshidx(length);
  for (int n=length-1; n>=1; n--) {
    int index=idx[n];
    threshidx[index] = length-n;
  }

  // Assign label to output according to size (descending)
  int xsize = labelim.xsize();
  int ysize = labelim.ysize();
  int zsize = labelim.zsize();
  for (int k=0; k<zsize; k++)
    for (int j=0; j<ysize; j++)
      for (int i=0; i<xsize; i++)
        outputimage[k*imgdim[0]*imgdim[1] + j*imgdim[0] + i] = threshidx[labelim.value(i, j, k)];

  output_nlabel[0] = length;

  return 0;

}


