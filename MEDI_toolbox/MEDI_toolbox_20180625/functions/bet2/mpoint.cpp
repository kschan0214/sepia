/*  Copyright (C) 1999-2004 University of Oxford  */

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

#include "mpoint.h"
#include "triangle.h"
#include "point.h"

namespace mesh{

Mpoint::~Mpoint(void){
}

Mpoint::Mpoint(double x, double y, double z, int counter, float val): _no(counter), _value(val) {
  _coord=Pt(x, y, z);
  _update_coord=Pt(0, 0, 0);
}

Mpoint::Mpoint(const Pt p, int counter,float val):_no(counter),_value(val) {
  _coord = p;
  _update_coord=Pt(0, 0, 0);
}


void Mpoint::translation(const double x,const double y,const double z)
{
  _coord+= Pt(x, y, z);
}
void Mpoint::rotation(const double r11, const double r12, const double r13, const double r21, const double r22, const double r23, const double r31, const double r32, const double r33,const double x, const double y, const double z) 
{
	Vec cen=_coord - Pt(x, y, z);
	
  _coord = Pt(x, y, z) + Vec( (cen.X*r11+cen.Y*r12+cen.Z*r13) , (cen.X*r21+cen.Y*r22+cen.Z*r23) , (cen.X*r31+cen.Y*r32+cen.Z*r33));
  
  
  
}
void Mpoint::rescale(const double t, const double x, const double y, const double z) 
{
  _coord = Pt(x, y, z) + t*(_coord - Pt(x, y, z));
}

void Mpoint::update() {
  _coord = _update_coord;
}

const Vec Mpoint::local_normal() const
{
  Vec v(0, 0, 0);
  for (list<Triangle*>::const_iterator i = _triangles.begin(); i!=_triangles.end(); i++)
    {
      v+=(*i)->normal();
    }
  v.normalize();
  return v;
}

const Pt Mpoint::medium_neighbours() const
{
  Pt resul(0, 0, 0);
  int counter=_neighbours.size();
  for (list<Mpoint*>::const_iterator i = _neighbours.begin(); i!=_neighbours.end(); i++)
    {
      resul+=(*i)->_coord;
    }
  resul=Pt(resul.X/counter, resul.Y/counter, resul.Z/counter);
  return resul;
}

const Vec Mpoint::difference_vector() const
{
  return medium_neighbours() - _coord;
}

const Vec Mpoint::orthogonal() const
{
  Vec n = local_normal();
  return n * (difference_vector()| n);
}

const Vec Mpoint::tangential() const
{
  return (difference_vector() - orthogonal());
}

const double Mpoint::medium_distance_of_neighbours() const
{
  double l = 0;
  for (list<Mpoint*>::const_iterator i=_neighbours.begin(); i!=_neighbours.end(); i++)
    {
      l+=(*(*i)-*this).norm();
    }
  l/=_neighbours.size();
  return l;
}
const Vec Mpoint::max_triangle() const
{  
  //returns a vector pointing from the vertex to the triangle centroid, scaled by the triangle area
 
  vector<float> Areas;
  int ind=0;
  Vec vA,temp;
  for (list<Triangle*>::const_iterator i=_triangles.begin(); i!=_triangles.end(); i++)
    {  
      temp=(*i)->area(this);
      Areas.push_back(temp.norm());
		//don't need to store in vector anymore
      if (Areas.back() >= Areas.at(ind)){
		ind=Areas.size()-1;
		vA = temp;
      }
    }

  return vA; 
}


const bool operator ==(const Mpoint &p2, const Mpoint &p1){
  return (fabs(p1.get_coord().X- p2.get_coord().X)<1e-8 && fabs(p1.get_coord().Y - p2.get_coord().Y)<1e-8 && fabs(p1.get_coord().Z - p2.get_coord().Z)<1e-8);
}

const Vec operator -(const Mpoint&p1, const Mpoint &p2){
  return Vec (p1.get_coord().X - p2.get_coord().X,p1.get_coord().Y - p2.get_coord().Y,p1.get_coord().Z - p2.get_coord().Z );
}

const Vec operator -(const Pt&p1, const Mpoint &p2){
  return Vec (p1.X - p2.get_coord().X,p1.Y - p2.get_coord().Y,p1.Z - p2.get_coord().Z );
}

const Vec operator -(const Mpoint&p1, const Pt &p2){
  return Vec (p1.get_coord().X - p2.X,p1.get_coord().Y - p2.Y,p1.get_coord().Z - p2.Z );
}

const bool operator <(const Mpoint &p1,const Mpoint &p2){
  bool result = false;
  for (list<Mpoint *>::const_iterator i= p1._neighbours.begin(); i!=p1._neighbours.end();i++){
    if (*(*i)==p2) result = true;
  } 
  return result;
}




}


