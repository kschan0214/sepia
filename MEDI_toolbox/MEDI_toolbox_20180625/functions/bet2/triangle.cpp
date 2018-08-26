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

#include "triangle.h"
#include "mpoint.h"

namespace mesh{

//this constructor also puts the connexions between the points. 
Triangle::Triangle(Mpoint* const p1, Mpoint* const p2, Mpoint* const p3,float val):_value(val) {
  oriented = false;
  _vertice[0]=p1;
  _vertice[1]=p2;
  _vertice[2]=p3;

  p1->_triangles.push_back(this);
  p2->_triangles.push_back(this);
  p3->_triangles.push_back(this);
  
  p1->_neighbours.remove(p2);
  p1->_neighbours.remove(p3);
  p2->_neighbours.remove(p3);
  p2->_neighbours.remove(p1);
  p3->_neighbours.remove(p1);
  p3->_neighbours.remove(p2);

  p1->_neighbours.push_back(p2);
  p1->_neighbours.push_back(p3);
  p2->_neighbours.push_back(p3);
  p2->_neighbours.push_back(p1);
  p3->_neighbours.push_back(p1);
  p3->_neighbours.push_back(p2);

}


//warning, you should remove neighbourhood relations between points by hand
Triangle::~Triangle() {
  _vertice[0]->_triangles.remove(this);
  _vertice[1]->_triangles.remove(this);
  _vertice[2]->_triangles.remove(this);
}


const Pt Triangle::centroid() const{
return Pt((_vertice[0]->get_coord().X +_vertice[1]->get_coord().X +_vertice[2]->get_coord().X)/3,
(_vertice[0]->get_coord().Y +_vertice[1]->get_coord().Y +_vertice[2]->get_coord().Y)/3,
(_vertice[0]->get_coord().Z +_vertice[1]->get_coord().Z +_vertice[2]->get_coord().Z)/3
);
}

const Vec Triangle::normal() const{
  Vec result = (_vertice[2]->get_coord() - _vertice[0]->get_coord()) * (_vertice[1]->get_coord() - _vertice[0]->get_coord());
  return result; 
}

const Vec Triangle::area(const Mpoint* const p) const{
  Vec v1,v2,vA;
  float Tarea;
 
  //calculate
  v1=*_vertice[1]-*_vertice[0];
  v2=*_vertice[2]-*_vertice[0];
  Tarea=0.5*((v1*v2).norm());
  //find appriopriate vector
  for (int i = 0; i<3; i++){
    if (p==_vertice[i]){
      vA=(this->centroid())-*_vertice[i];
    }
  }
  vA=vA/vA.norm()*Tarea;

  return vA; 
}



Mpoint * const Triangle::get_vertice(const int i) const
{return _vertice[i];}

void Triangle::swap() {
  Mpoint * p = _vertice[1];
  _vertice[1] = _vertice[2];
  _vertice[2] = p;
}

//check if two triangles are adjacents
//0 if not
//1 if yes and good orientation
//2 if yes and bad orientation
const int Triangle::operator <(const Triangle * const t) const{
  int c = 0;
  int a11=-1, a12=-1, a21=-1, a22=-1;
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      if (_vertice[i]==t->_vertice[j])
	{
	  if (a11 == -1) {a11=i; a21=j;}
	  else {a12=i; a22=j;};
	  c++;//cout<<i<<"++"<<j<<endl;
	};
  if (c == 2) {if ((a12-a11 + a22-a21) % 3 == 0) return 1;
  else return 2;}
  else return 0;
}

const bool Triangle::operator ==(const Triangle & t) const
{
  return ((*get_vertice(0)) == *(t.get_vertice(0)) && (*get_vertice(1)) == *(t.get_vertice(1)) && (*get_vertice(2)) == *(t.get_vertice(2)));
}

const bool Triangle::intersect(const Triangle & t) const
{
  bool result = false;
  Vec normal = (this->get_vertice(0)->get_coord() - this->get_vertice(1)->get_coord()) * (this->get_vertice(0)->get_coord() - this->get_vertice(2)->get_coord());

  if ((normal|(t.get_vertice(0)->get_coord() - this->get_vertice(0)->get_coord())) * (normal|(t.get_vertice(1)->get_coord() - this->get_vertice(0)->get_coord())) < 0)
    {
      //possible intersection -> make the full test
      //test from *this
      for (int i = 0; i < 3; i++)
	if ((normal|(t.get_vertice(i)->get_coord() - this->get_vertice(0)->get_coord())) * (normal|(t.get_vertice((i+1)%3)->get_coord() - this->get_vertice(0)->get_coord())) < 0)
	  {
	    Vec v1 = this->get_vertice(1)->get_coord() - this->get_vertice(0)->get_coord(); 
	    Vec v2 = this->get_vertice(2)->get_coord() - this->get_vertice(0)->get_coord(); 
	    Vec v3 = this->get_vertice(2)->get_coord() - this->get_vertice(1)->get_coord(); 
	    Vec v = v1 * v2;
	    
	    Vec p1 = t.get_vertice(i)->get_coord() - this->get_vertice(0)->get_coord();
	    Vec d1 = t.get_vertice((i+1)%3)->get_coord() - t.get_vertice(i)->get_coord();
	    double denom = (d1.X * v.X + d1.Y * v.Y + d1.Z * v.Z);
	    if (denom != 0)
	      {
		double lambda1 = - (p1.X * v.X + p1.Y * v.Y + p1.Z * v.Z)/denom;
		Vec proj1 = p1 + (d1 * lambda1);
		
		//checks if proj is inside the triangle ...
		bool inside = false;
		Vec n1 = v1 * proj1;
		Vec n2 = proj1 * v2;
		Vec n3 = v3 * (proj1 + (v1 * -1));
		if (((n1 | n3) > 0 & (n2 | n3) > 0 & (n1 | n2) > 0) | ((n1 | n3) < 0 & (n2 | n3) < 0 & (n1 | n2) < 0) )
		  inside = true;
		
		result = result | inside;
	      }
	  }

      //test from t
      
      Vec normalt = (t.get_vertice(0)->get_coord() - t.get_vertice(1)->get_coord()) * (t.get_vertice(0)->get_coord() - t.get_vertice(2)->get_coord());
      for (int i = 0; i < 3; i++)
	if ((normalt|(this->get_vertice(i)->get_coord() - t.get_vertice(0)->get_coord())) * (normalt|(this->get_vertice((i+1)%3)->get_coord() - t.get_vertice(0)->get_coord())) < 0)
	  {
	    Vec v1 = t.get_vertice(1)->get_coord() - t.get_vertice(0)->get_coord(); 
	    Vec v2 = t.get_vertice(2)->get_coord() - t.get_vertice(0)->get_coord(); 
	    Vec v3 = t.get_vertice(2)->get_coord() - t.get_vertice(1)->get_coord(); 
	    Vec v = v1 * v2;
	    
	    Vec p1 = this->get_vertice(i)->get_coord() - t.get_vertice(0)->get_coord();
	    Vec d1 = this->get_vertice((i+1)%3)->get_coord() - this->get_vertice(i)->get_coord();

	    double denom = (d1.X * v.X + d1.Y * v.Y + d1.Z * v.Z);
	    if (denom != 0)
	      {
		double lambda1 = - (p1.X * v.X + p1.Y * v.Y + p1.Z * v.Z)/denom;
		Vec proj1 = p1 + (d1 * lambda1);
		
		//checks if proj is inside the triangle ...
		bool inside = false;
		Vec n1 = v1 * proj1;
		Vec n2 = proj1 * v2;
		Vec n3 = v3 * (proj1 + (v1 * -1));
		if (((n1 | n3) > 0 & (n2 | n3) > 0 & (n1 | n2) > 0) | ((n1 | n3) < 0 & (n2 | n3) < 0 & (n1 | n2) < 0) )
		  inside = true;
		
		result = result | inside;
	      }
	  }



    }
  else if ((normal|(t.get_vertice(0)->get_coord() - this->get_vertice(0)->get_coord())) * (normal|(t.get_vertice(2)->get_coord() - this->get_vertice(0)->get_coord())) < 0)
    {
      //possible intersection -> make the full test
      //test from *this
      for (int i = 0; i < 3; i++)
	if ((normal|(t.get_vertice(i)->get_coord() - this->get_vertice(0)->get_coord())) * (normal|(t.get_vertice((i+1)%3)->get_coord() - this->get_vertice(0)->get_coord())) < 0)
	  {
	    Vec v1 = this->get_vertice(1)->get_coord() - this->get_vertice(0)->get_coord(); 
	    Vec v2 = this->get_vertice(2)->get_coord() - this->get_vertice(0)->get_coord(); 
	    Vec v3 = this->get_vertice(2)->get_coord() - this->get_vertice(1)->get_coord(); 
	    Vec v = v1 * v2;
	    
	    Vec p1 = t.get_vertice(i)->get_coord() - this->get_vertice(0)->get_coord();
	    Vec d1 = t.get_vertice((i+1)%3)->get_coord() - t.get_vertice(i)->get_coord();
	    double denom = (d1.X * v.X + d1.Y * v.Y + d1.Z * v.Z);
	    if (denom != 0)
	      {
		double lambda1 = - (p1.X * v.X + p1.Y * v.Y + p1.Z * v.Z)/denom;
		Vec proj1 = p1 + (d1 * lambda1);
		//checks if proj is inside the triangle ...
		bool inside = false;
		Vec n1 = v1 * proj1;
		Vec n2 = proj1 * v2;
		Vec n3 = v3 * (proj1 + (v1 * -1));
		if (((n1 | n3) > 0 & (n2 | n3) > 0 & (n1 | n2) > 0) | ((n1 | n3) < 0 & (n2 | n3) < 0 & (n1 | n2) < 0) )
		  {
		    inside = true;
		  }
		result = result | inside;
	      }
	  }
      
      //test from t
      Vec normalt = (t.get_vertice(0)->get_coord() - t.get_vertice(1)->get_coord()) * (t.get_vertice(0)->get_coord() - t.get_vertice(2)->get_coord());
      for (int i = 0; i < 3; i++)
	if ((normalt|(this->get_vertice(i)->get_coord() - t.get_vertice(0)->get_coord())) * (normalt|(this->get_vertice((i+1)%3)->get_coord() - t.get_vertice(0)->get_coord())) < 0)
	  {
	    Vec v1 = t.get_vertice(1)->get_coord() - t.get_vertice(0)->get_coord(); 
	    Vec v2 = t.get_vertice(2)->get_coord() - t.get_vertice(0)->get_coord(); 
	    Vec v3 = t.get_vertice(2)->get_coord() - t.get_vertice(1)->get_coord(); 
	    Vec v = v1 * v2;
	    
	    Vec p1 = this->get_vertice(i)->get_coord() - t.get_vertice(0)->get_coord();
	    Vec d1 = this->get_vertice((i+1)%3)->get_coord() - this->get_vertice(i)->get_coord();
	    double denom = (d1.X * v.X + d1.Y * v.Y + d1.Z * v.Z);
	    if (denom != 0)
	      {
		double lambda1 = - (p1.X * v.X + p1.Y * v.Y + p1.Z * v.Z)/denom;
		Vec proj1 = p1 + (d1 * lambda1);
		
		//checks if proj is inside the triangle ...
		bool inside = false;
		Vec n1 = v1 * proj1;
		Vec n2 = proj1 * v2;
		Vec n3 = v3 * (proj1 + (v1 * -1));
		if (((n1 | n3) > 0 & (n2 | n3) > 0 & (n1 | n2) > 0) | ((n1 | n3) < 0 & (n2 | n3) < 0 & (n1 | n2) < 0) )
		  {
		    inside = true;
		  }
		
		result = result | inside;
	      }
	  }




    }
  else {return (false);}
  return result;



}


}



