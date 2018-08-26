
/*  BET - Brain Extraction Tool

    BETv1 Steve Smith
    BETv2 Mickael Pechaud, Mark Jenkinson, Steve Smith
    FMRIB Image Analysis Group

    Copyright (C) 1999-2006 University of Oxford  */

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

#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <algorithm>

#include "utils/options.h"
#include "newimage/newimageall.h"
#include "meshclass/meshclass.h"

#include "bet2.h"
#include "mex.h"

using namespace std;
using namespace NEWIMAGE;
using namespace Utilities;
using namespace mesh;

void noMoreMemory()
{
  cerr<<"Unable to satisfy request for memory"<<endl;
  abort();
}

struct bet_parameters
{
  double min, max, t98, t2, t, tm, radius;
  Pt cog;
};

const double normal_max_update_fraction = .5;
const double lambda_fit = .1;


vector<float> empty_vector(0, 0);

string title="BET (Brain Extraction Tool) v2.1 - FMRIB Analysis Group, Oxford";
string examples="bet2 <input_fileroot> <output_fileroot> [options]";

Option<bool> verbose(string("-v,--verbose"), false, 
		     string("switch on diagnostic messages"), 
		     false, no_argument);
Option<bool> generate_mesh(string("-e,--mesh"), false, 
		     string("generates brain surface as mesh in vtk format"), 
		     false, no_argument);
Option<bool> help(string("-h,--help"), false, 
		     string("displays this help, then exits"), 
		     false, no_argument);
Option<bool> outline(string("-o,--outline"), false, 
		     string("generate brain surface outline overlaid onto original image"), 
		     false, no_argument);
Option<bool> skull(string("-s,--skull"), false, 
		   string("generate approximate skull image"),
		   false, no_argument);
Option<bool> mask(string("-m,--mask"), false, 
		  string("generate binary brain mask"), 
		  false, no_argument);
Option<bool> no_output(string("-n,--nooutput"), false, 
		       string("don't generate segmented brain image output"), 
		       false, no_argument);
Option<bool> apply_thresholding(string("-t,--threshold"), false, 
				string("-apply thresholding to segmented brain image and mask"), 
				false, no_argument);
Option<float> fractional_threshold(string("-f"), 0.5,
				    string("~<f>\t\tfractional intensity threshold (0->1); default=0.5; smaller values give larger brain outline estimates"), false, requires_argument);

Option<float> gradient_threshold(string("-g"), 0.0,
				 string("~<g>\t\tvertical gradient in fractional intensity threshold (-1->1); default=0; positive values give larger brain outline at bottom, smaller at top"), false, requires_argument);

Option<float> radiusarg(string("-r,--radius"), 0.0,
		     string("~<r>\thead radius (mm not voxels); initial surface sphere is set to half of this"), false, requires_argument);

Option<float> centerarg(string("-c"), 0.0,
				 string("~<x y z>\tcentre-of-gravity (voxels not mm) of initial mesh surface."), false, requires_3_arguments);


void draw_segment(volume<short>& image, const Pt& p1, const Pt& p2)
{
  double xdim = (double) image.xdim();
  double ydim = (double) image.ydim();
  double zdim = (double) image.zdim();
  double mininc = min(xdim,min(ydim,zdim)) * .5;

  Vec n = p1 - p2;
  double d = n.norm();
  n.normalize();

  for (double i=0; i<=d; i+=mininc)
    {
      Pt p = p2 + i* n;
      image((int) floor((p.X)/xdim +.5),(int) floor((p.Y)/ydim +.5),(int) floor((p.Z)/zdim +.5)) = 0;
    }
}


volume<short> draw_mesh(const volume<short>& image, const Mesh &m)
{
  double xdim = (double) image.xdim();
  double ydim = (double) image.ydim();
  double zdim = (double) image.zdim();
  double mininc = min(xdim,min(ydim,zdim)) * .5;
  volume<short> res = image;
  for (list<Triangle*>::const_iterator i = m._triangles.begin(); i!=m._triangles.end(); i++)
    {
      Vec n = (*(*i)->get_vertice(0) - *(*i)->get_vertice(1));
      double d = n.norm();
      n.normalize();

      for (double j=0; j<=d; j+=mininc)
	{
	  Pt p = (*i)->get_vertice(1)->get_coord()  + j* n;
	  draw_segment(res, p, (*i)->get_vertice(2)->get_coord());
	} 
    }
  return res;
}

void draw_segment_bis(volume<float>& image, const Pt& p1, const Pt& p2, double max)
{
  double xdim = (double) image.xdim();
  double ydim = (double) image.ydim();
  double zdim = (double) image.zdim();
  double mininc = min(xdim,min(ydim,zdim)) * .5;

  Vec n = p1 - p2;
  double d = n.norm();
  n.normalize();

  for (double i=0; i<=d; i+=mininc)
    {
      Pt p = p2 + i* n;
      image((int) ((p.X)/xdim +.5),(int) ((p.Y)/ydim +.5),(int) ((p.Z)/zdim +.5)) = max;
    }
}


volume<float> draw_mesh_bis(const volume<float>& image, const Mesh &m)
{
  double xdim = (double) image.xdim();
  double ydim = (double) image.ydim();
  double zdim = (double) image.zdim();
  double mininc = min(xdim,min(ydim,zdim)) * .5;
  volume<float> res = image;
  double max = image.max();
  for (list<Triangle*>::const_iterator i = m._triangles.begin(); i!=m._triangles.end(); i++)
    {
      Vec n = (*(*i)->get_vertice(0) - *(*i)->get_vertice(1));
      double d = n.norm();
      n.normalize();

      for (double j=0; j<=d; j+=mininc)
	{
	  Pt p = (*i)->get_vertice(1)->get_coord()  + j* n;
	  draw_segment_bis(res, p, (*i)->get_vertice(2)->get_coord(), max);
	} 
    }
  return res;
}



//calculates the mask from the mesh by spreading an initial point outside the mesh, and stopping it when the mesh is reached.
volume<short> make_mask_from_mesh(const volume<float> & image, const Mesh& m)
{
  //  cout<<"make_mask_from_mesh begins"<<endl;

  double xdim = (double) image.xdim();
  double ydim = (double) image.ydim();
  double zdim = (double) image.zdim();

  volume<short> mask;
  copyconvert(image,mask);
  
  int xsize = mask.xsize();
  int ysize = mask.ysize();
  int zsize = mask.zsize();
  
  mask = 1;
  
  mask = draw_mesh(mask, m);

  vector<Pt> current;
  current.clear();
  Pt c(0., 0., 0.);
  for (vector<Mpoint *>::const_iterator it=m._points.begin(); it!=m._points.end(); it++)
    c+=(*it)->get_coord();

  c*=(1./m._points.size());
  c.X/=xdim; c.Y/=ydim; c.Z/=zdim;

  current.push_back(c);

  while (!current.empty())
    {
      Pt pc = current.back();
      int x, y, z;
      x=(int) pc.X; y=(int) pc.Y; z=(int) pc.Z;
      //current.remove(pc);
      current.pop_back();
      mask.value(x, y, z) = 0;
      if (0<=x-1 && mask.value(x-1, y, z)==1) current.push_back(Pt(x-1, y, z));
      if (0<=y-1 && mask.value(x, y-1, z)==1) current.push_back(Pt(x, y-1, z));
      if (0<=z-1 && mask.value(x, y, z-1)==1) current.push_back(Pt(x, y, z-1));
      if (xsize>x+1 && mask.value(x+1, y, z)==1) current.push_back(Pt(x+1, y, z));
      if (ysize>y+1 && mask.value(x, y+1, z)==1) current.push_back(Pt(x, y+1, z));
      if (zsize>z+1 && mask.value(x, y, z+1)==1) current.push_back(Pt(x, y, z+1)); 
    }

  //  cout<<"make_mask_from_mesh ends"<<endl;
  return mask;
}


volume<float> inside_mesh(const volume<float> & image, const Mesh& m)
{
  volume<float> res = image;
  int xsize = image.xsize();
  int ysize = image.ysize();
  int zsize = image.zsize();
  volume<short> inside = make_mask_from_mesh(image, m);
  for (int k=0; k<zsize; k++)
    for (int j=0; j<ysize; j++)
      for (int i=0; i<xsize; i++)
	res.value(i, j, k) = (1-inside.value(i, j, k)) * image.value(i, j, k);
  return res;
}



bet_parameters adjust_initial_mesh(const volume<float> & image, Mesh& m, const double & rad = 0., const double xpara=0.,  const double ypara=0.,  const double zpara=0.)
{
  bet_parameters bp;
  double xdim = image.xdim();
  double ydim = image.ydim();
  double zdim = image.zdim();
  double t2, t98, t;

  //computing t2 && t98
  //  cout<<"computing robust min && max begins"<<endl;

  bp.min = image.min();
  bp.max = image.max();

  t2 = image.robustmin();
  t98 = image.robustmax();
  //t2=32.;
  //t98=16121.;
  
  //  cout<<"computing robust min && max ends"<<endl;
  
  t = t2 + .1*(t98 - t2);
  bp.t98 = t98;
  bp.t2 = t2;
  bp.t = t;
  //  cout<<"t2 "<<t2<<" t98 "<<t98<<" t "<<t<<endl;
  
  //  cout<<"computing center && radius begins"<<endl;
  
  //finds the COG
  Pt center(0, 0, 0);
  double counter = 0;
  if ((xpara == 0.) & (ypara==0.) & (zpara==0.))
    {
      double tmp = t - t2;
      for (int k=0; k<image.zsize(); k++)
	for (int j=0; j<image.ysize(); j++)
	  for (int i=0; i<image.xsize(); i++)
	    {
	      double c = image(i, j, k ) - t2;
	      if (c > tmp)
		{
		  c = min(c, t98 - t2);   
		  counter+=c;
		  center +=  Pt(c*i*xdim, c*j*ydim, c*k*zdim);
		}
	    }
      center=Pt(center.X/counter, center.Y/counter, center.Z/counter);
      //cout<<counter<<endl;
      //  cout<<"cog "<<center.X<<" "<<center.Y<<" "<<center.Z<<endl;
    }
  else center=Pt(xpara*xdim, ypara*ydim, zpara*zdim);
  
  bp.cog = center;

  if (rad == 0.)
    {
      double radius=0;
      counter=0;
      double scale=xdim*ydim*zdim;
      for (int k=0; k<image.zsize(); k++)
	for (int j=0; j<image.ysize(); j++)
	  for (int i=0; i<image.xsize(); i++)
	    {
	      double c = image(i, j, k);
	      if (c > t)
		{
		  counter+=1;
		}
	    }
      radius = pow (.75 * counter*scale/M_PI, 1.0/3.0);
      //      cout<<radius<<endl;
      bp.radius = radius;
    } 
  else (bp.radius = rad);

  m.translation(center.X, center.Y, center.Z);
  m.rescale (bp.radius/2, center);

  //  cout<<"computing center && radius ends"<<endl;

  //computing tm
  //  cout<<"computing tm begins"<<endl;
  vector<double> vm;
  for (int k=0; k<image.zsize(); k++)
    for (int j=0; j<image.ysize(); j++)
      for (int i=0; i<image.xsize(); i++)
	{
	  double d = image.value(i, j, k);
	  Pt p(i*xdim, j*ydim, k*zdim);
	  if (d > t2 && d < t98 && ((p - center)|(p - center)) < bp.radius * bp.radius)
	    vm.push_back(d);
	}

  int med = (int) floor(vm.size()/2.);
  //  cout<<"before sort"<<endl;
  nth_element(vm.begin(), vm.begin() + med - 1, vm.end());
  //partial_sort(vm.begin(), vm.begin() + med + 1, vm.end());
  //double tm = vm[med];
  double tm=(*max_element(vm.begin(), vm.begin() + med));
  //  cout<<"tm "<<tm<<endl;
  bp.tm = tm;
  //  cout<<"computing tm ends"<<endl;
  
  return bp;
}



double step_of_computation(const volume<float> & image, Mesh & m, const double bet_main_parameter, const int pass, const double increase_smoothing, const int iteration_number, double & l, const double t2, const double tm, const double t, const double E,const double F, const double zcog, const double radius, const double local_th=0., const int d1=7, const int d2=3){
  double xdim = image.xdim();
  double ydim = image.ydim();
  double zdim = image.zdim();
  double dscale = Min(Min(Min(xdim,ydim),zdim),1.0);
  //cout << xdim << " " << ydim << " " << zdim << " " << dscale << endl;
  
  if (iteration_number==50 || iteration_number%100 == 0 )
    {
      l = 0;
      int counter = 0;
      for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ )
	{
	  counter++;
	  l += (*i)->medium_distance_of_neighbours();
	}
      l/=counter;
    }
  
  for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++)
    {
      Vec sn, st, u1, u2, u3, u;
      double f2, f3=0;
      
      Vec n = (*i)->local_normal();
      Vec dv = (*i)->difference_vector();
      
      double tmp = dv|n;
      sn = n * tmp;
      st = dv - sn;
      
      u1 = st*.5;
      
      double rinv = (2 * fabs(sn|n))/(l*l);
      
      f2 = (1+tanh(F*(rinv - E)))*0.5;
      if (pass > 0)
	if (tmp > 0) {
	  f2*=increase_smoothing;
	  f2 = Min(f2, 1.);
	}
      
      u2 = f2 * sn;
      
      //main term of bet
      {
	double local_t = bet_main_parameter;
	if (local_th != 0.0)
	  {
	    local_t = Min(1., Max(0., bet_main_parameter + local_th*((*i)->get_coord().Z - zcog)/radius));
	  }
	
	double Imin = tm;
	double Imax = t;
	
	Pt p = (*i)->get_coord() + (-1)*n;
	double iv = p.X/xdim + .5, jv = p.Y/ydim +.5, kv = p.Z/zdim +.5; 
	if (image.in_bounds((int)iv,(int) jv,(int) kv))
	  {	
	    double im=image.value((int)iv,(int)jv,(int)kv);
	    Imin = Min(Imin, im);
	    Imax = Max(Imax,im);
	    
	    double nxv=n.X/xdim, nyv=n.Y/ydim, nzv=n.Z/zdim;
	    int i2=(int)(iv-(d1-1)*nxv), j2 =(int) (jv-(d1-1)*nyv), k2 =(int)(kv-(d1-1)*nzv); 
	    nxv*=dscale; nyv*=dscale; nzv*=dscale;
	    if (image.in_bounds(i2, j2, k2))
	      {	
		im=image.value(i2,j2,k2);
		Imin = Min(Imin, im);
		
		for (double gi=2.0; gi<d1; gi+=dscale)
		  {
		    //cout << gi << " " << endl;
		    // the following is a quick calc of Pt p = (*i)->get_coord() + (-gi)*n;
		    iv-=nxv; jv-=nyv; kv-=nzv;
		    im = image.value((int) (iv), (int) (jv), (int) (kv));
		    Imin = Min(Imin, im);
		
		    if (gi<d2)
		      Imax = Max(Imax,im);
		  }
		
		Imin = Max (t2, Imin);
		Imax = Min (tm, Imax);	
		
		const double tl = (Imax - t2) * local_t + t2;
		
		if (Imax - t2 > 0)
		  f3=2*(Imin - tl)/(Imax - t2);
		else f3=(Imin - tl)*2;
	      }
	  }
	
      }
      
      f3 *= (normal_max_update_fraction * lambda_fit * l);
      
      u3 = f3 * n;
      
      u = u1 + u2 + u3;
            
      //cout<<"l "<<l<<"u1 "<<((u1*n).norm())<<"u2 "<<(u2|n)<<"u3 "<<(u3|n)<<endl;
      
      (*i)->_update_coord = (*i)->get_coord() + u;
    }

  m.update();
  
  return (0); 
}



volume<float> find_skull (volume<float> & image, const Mesh & m, const double t2, double t, double t98)
{ 
  const double skull_search = 30;
  const double skull_start = -3;
  
  volume<float> result = image;
  result=0;

  volume<short> volmesh;
  copyconvert(image,volmesh);  
  int xsize = volmesh.xsize();
  int ysize = volmesh.ysize();
  int zsize = volmesh.zsize();  
  double xdim = volmesh.xdim();
  double ydim = volmesh.ydim();
  double zdim = volmesh.zdim();  
  double scale = Min(xdim, Min(ydim, zdim)); 
  volmesh = 1;
  volmesh = draw_mesh(volmesh, m);
  
  image.setinterpolationmethod(trilinear);

  for (vector<Mpoint*>::const_iterator i = m._points.begin(); i != m._points.end(); i++)
    {
      double max_neighbour = 0;
      const Vec normal = (*i)->local_normal();
      const Vec n = Vec(normal.X/xdim, normal.Y/ydim, normal.Z/zdim);      

      for (list<Mpoint*>::const_iterator nei = (*i)->_neighbours.begin(); nei != (*i)->_neighbours.end(); nei++)
	max_neighbour = Max(((**i) - (**nei)).norm(), max_neighbour); 

      max_neighbour = ceil((max_neighbour)/2);

      const Pt mpoint((*i)->get_coord().X/xdim,(*i)->get_coord().Y/ydim,(*i)->get_coord().Z/zdim);
      for (int ck = (int)floor(mpoint.Z - max_neighbour/zdim); ck <= (int)floor(mpoint.Z + max_neighbour/zdim); ck++)
	for (int cj = (int)floor(mpoint.Y - max_neighbour/ydim); cj <= (int)floor(mpoint.Y + max_neighbour/ydim); cj++)
	  for (int ci = (int)floor(mpoint.X - max_neighbour/xdim); ci <= (int)floor(mpoint.X + max_neighbour/xdim); ci++)
	    {
	      bool compute = false;
	      const Pt point(ci, cj, ck);
	      const Pt realpoint(ci*xdim, cj*ydim, ck*zdim);
	      if (volmesh(ci, cj, ck) == 0) 
		{
		  double mindist = 10000;
		  for (list<Mpoint*>::const_iterator nei = (*i)->_neighbours.begin(); nei != (*i)->_neighbours.end(); nei++)
		    mindist = Min(((realpoint) - (**nei)).norm(), mindist); 
		  if (mindist >= ((realpoint) - (**i)).norm()) compute = true;
		}
	    

	      if (compute)
		{
		  double maxval = t;
		  double minval = image.interpolate(point.X, point.Y, point.Z);
		  double d_max = 0;
		  
		  for (double d=0; d<skull_search; d+=scale*.5)
		    {
		      Pt current = point + d * n;
		      double val = image.interpolate(current.X, current.Y, current.Z);
		      if (val>maxval)
			{
			  maxval=val;
			  d_max=d;
			}
		      
		      if (val<minval)
			minval=val;
		    }
		  
		  if (maxval > t)
		    {
		      double d_min=skull_start;
		      double maxJ =-1000000;
		      double lastJ=-2000000;
		      for(double d=skull_start; d<d_max; d+=scale*0.5)
			{
			  Pt current = point + d * n;
			  if (current.X >= 0 && current.Y >= 0 && current.Z >= 0 && current.X<xsize && current.Y<ysize && current.Z<zsize)
			    {
			      double tmpf = d/30 - image.interpolate(current.X, current.Y, current.Z) / (t98 - t2);
			      if (tmpf > maxJ)
				{
				  maxJ=tmpf;
				  d_min = d;
				}
			      lastJ=tmpf;
			    }
			}
		      double maxgrad = 0;
		      double d_skull;
		      Pt current2 = point + d_min * n;
		      if (current2.X >= 0 && current2.Y >= 0 && current2.Z >= 0 && current2.X<xsize && current2.Y<ysize && current2.Z<zsize)
			{
			  double val2 = image.interpolate(current2.X, current2.Y, current2.Z);
			  for(double d=d_min + scale; d<d_max; d+=0.5*scale)
			    {
			      Pt current = point + d * n;
			      if (current.X >= 0 && current.Y >= 0 && current.Z >= 0 && current.X<xsize && current.Y<ysize && current.Z<zsize)
				{
				  double val = image.interpolate(current.X, current.Y, current.Z);
				  double grad = val - val2;
				  val2 = val;
				  if (grad > 0)
				    {
				      if (grad > maxgrad)
					{
					  maxgrad=grad;
					  d_skull=d;
					}
				      else d = d_max;
				    }
				}
			    }
			}
		      if (maxgrad > 0)
			{
			  Pt current3 = point + d_skull * n;
			  if (current3.X >= 0 && current3.Y >= 0 && current3.Z >= 0 && current3.X<xsize && current3.Y<ysize && current3.Z<zsize)
			    result ((int)current3.X, (int)current3.Y, (int)current3.Z) = 100/*max*/;
			}
		    }
		}
	    }
    }
  return result;
}

template <typename T>
int bet_impl(float *inputimage, T outputmask, const int* imgdim, const float *voxel_size)
{

  //set a memory hanlder that displays an error message
  set_new_handler(noMoreMemory);


  //the real program

  volume<float> testvol(imgdim[0], imgdim[1], imgdim[2], inputimage, 0);
  testvol.setdims(voxel_size[0],voxel_size[1], voxel_size[2]); 
  
//  if (read_volume(testvol,in.c_str())<0)  return -1;

  double xarg = 0, yarg = 0, zarg = 0;
  if (centerarg.set())
    {
      /*
      if (centerarg.value().size()!=3)
	{
	  cerr<<"three parameters expected for center option !"<<endl;
	  cerr<<"please check there is no space after commas."<<endl;
	  exit (-1);
	}
      else
      */
	{
	  xarg = (double) centerarg.value(0);
	  yarg = (double) centerarg.value(1);
	  zarg = (double) centerarg.value(2);
	  ColumnVector v(4);
	  v << xarg << yarg << zarg << 1.0;
	  v = testvol.niftivox2newimagevox_mat() * v;
	  xarg = v(1);  yarg = v(2);  zarg = v(3);
	}
    }
  
  const double bet_main_parameter = pow(fractional_threshold.value(), .275);
  

  // 2D kludge (worked for bet, but not here in bet2, hohum)
  if (testvol.xsize()*testvol.xdim()<20) testvol.setxdim(200);
  if (testvol.ysize()*testvol.ydim()<20) testvol.setydim(200);
  if (testvol.zsize()*testvol.zdim()<20) testvol.setzdim(200);
  
  Mesh m;
  make_mesh_from_icosa(5, m); 
  
  
  bet_parameters bp = adjust_initial_mesh(testvol, m, radiusarg.value(), xarg, yarg, zarg);
  
  if (verbose.value())
    {
      cout<<"min "<<bp.min<<" thresh2 "<<bp.t2<<" thresh "<<bp.t<<" thresh98 "<<bp.t98<<" max "<<bp.max<<endl;
      cout<<"c-of-g "<<bp.cog.X<<" "<<bp.cog.Y<<" "<<bp.cog.Z<<" mm"<<endl;
      cout<<"radius "<<bp.radius<<" mm"<<endl;
      cout<<"median within-brain intensity "<<bp.tm<<endl;
    }
  

  Mesh moriginal=m;
  
  const double rmin=3.33;
  const double rmax=10;
  const double E = (1/rmin + 1/rmax)/2.;
  const double F = 6./(1/rmin - 1/rmax);
  const int nb_iter = 1000;
  const double self_intersection_threshold = 4000;
  
  double l = 0;
  for (int i=0; i<nb_iter; i++)
    {
      step_of_computation(testvol, m, bet_main_parameter, 0, 0, i, l, bp.t2, bp.tm, bp.t, E, F, bp.cog.Z, bp.radius, gradient_threshold.value());
    }
  
  double tmp = m.self_intersection(moriginal);
  if (verbose.value() && !generate_mesh.value())
    cout<<"self-intersection total "<<tmp<<" (threshold=4000.0) "<<endl;
  
  bool self_intersection;
  if (!generate_mesh.value()) self_intersection = (tmp > self_intersection_threshold);  
  else (self_intersection = m.real_self_intersection());
  int pass = 0;  

  if (verbose.value() && generate_mesh.value() && self_intersection)
    cout<<"the mesh is self-intersecting "<<endl;


  //self-intersection treatment
  while (self_intersection)
    {
      if (self_intersection && verbose.value()) {cout<<"thus will rerun with higher smoothness constraint"<<endl;};
      m = moriginal;
      l = 0;
      pass++;
      for (int i=0; i<nb_iter; i++)
	{
	  double incfactor = pow (10.0,(double) pass + 1);
	  if (i > .75 * (double)nb_iter)
	    incfactor = 4.*(1. - i/(double)nb_iter) * (incfactor - 1.) + 1.;
	  step_of_computation(testvol, m, bet_main_parameter, pass, incfactor, i, l, bp.t2, bp.tm, bp.t, E, F, bp.cog.Z, bp.radius, gradient_threshold.value());
	}
      double tmp = m.self_intersection(moriginal);
  
      self_intersection = (tmp > self_intersection_threshold);
      if (!generate_mesh.value()) self_intersection = (tmp > self_intersection_threshold);  
      else (self_intersection = m.real_self_intersection());
      
      if (verbose.value() && !generate_mesh.value())
	cout<<"self-intersection total "<<tmp<<" (threshold=4000.0) "<<endl;
      
      if (verbose.value() && generate_mesh.value() && self_intersection)
	cout<<"the mesh is self-intersecting"<<endl;
	
      if (pass==10) // give up
	self_intersection=0;
    }
  

  //display
  volume<short> brainmask = make_mask_from_mesh(testvol, m);
  
	{	
    	int xsize = testvol.xsize();
      	int ysize = testvol.ysize();
      	int zsize = testvol.zsize();
      	for (int k=0; k<zsize; k++)
			for (int j=0; j<ysize; j++)
	  			for (int i=0; i<xsize; i++)
	    			outputmask[k*imgdim[0]*imgdim[1] + j*imgdim[0] + i] = 1 - brainmask.value(i, j, k);
    }
  

  return 0;
  
}

int bet(float *inputimage, int* outputmask, const int* imgdim, const float *voxel_size)
{
  return bet_impl<int*>(inputimage, outputmask, imgdim, voxel_size);
}

// ------------------- interface ----------------------------

void mexFunction( int nlhs, mxArray      *plhs[],
		          int nrhs, const mxArray *prhs[]) {
    if (nrhs != 3) {
    mexErrMsgTxt("MEX requires 3 inputs."); }
    
    //int N1, N2, N3, peel, depth;
    //int *mask, *m_size, n;
    //float* fT, *v_size, tol;
    int m_size[3];
    float v_size[3];
        
    //mexPrintf("%g %d %d\n",tol,depth,N1);
    //mexEvalString("drawnow;");

    int i;
    double *p,*pp;

    p = mxGetPr(prhs[1]);
    pp = mxGetPr(prhs[2]);
    for (i=0;i<3;++i){
        m_size[i] = p[i];
        v_size[i] = pp[i];
    }

mexPrintf("%d %d %d\n",m_size[0],m_size[1],m_size[2]);
mexPrintf("%g %g %g\n",v_size[0],v_size[1],v_size[2]);

    int n = mxGetNumberOfElements(prhs[0]);
mexPrintf("%d \n",n);
mexEvalString("drawnow;");

    float *fM = new float [n]();
    int *Mask = new int [n]();
    p = mxGetPr(prhs[0]);
    for (i=0;i<n;++i){
        fM[i] = p[i];
    }

    bet(fM, Mask, m_size, v_size);

    mxArray *out_m;
    out_m = plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
    double *c;
    c = mxGetPr(out_m);
    for (int i=0;i<n;++i)
        c[i] = Mask[i];

    // free memory
    delete fM;
    delete Mask;

}









