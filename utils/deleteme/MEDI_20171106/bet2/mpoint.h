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

#ifndef _mpoint
#define _mpoint

#include <list>
#include <iostream>

using namespace std;

#include "point.h"

namespace mesh{

class Triangle;

class Mpoint {
 public:
  ~Mpoint(void);
  Mpoint(double x, double y, double z, int counter,float val=0);
  Mpoint(const Pt p, int counter,float val=0);

  Mpoint (Mpoint & m); //prevents from using copy constructor
  Mpoint operator=(Mpoint & m); //prevents from using affectation operator

  void translation(const double x,const double y,const double z);
  void rescale(const double t, const double x, const double y, const double z);
  void rotation(const double r11, const double r12, const double r13, const double r21, const double r22, const double r23, const double r31, const double r32, const double r33,const double x, const double y, const double z);

  void update();
  
  Pt _update_coord;
  list<Triangle*> _triangles;

  const Vec local_normal() const;  
  const Pt medium_neighbours() const;
  const Vec difference_vector() const;
  const Vec orthogonal() const;
  const Vec tangential() const;
  const Vec max_triangle() const;	
  const double medium_distance_of_neighbours() const;

  const Pt get_coord() const {return _coord;};
  const int get_no() const {return _no;};
  const float get_value() const {return _value;};
  void set_value(const float val){_value=val;};
  list<Mpoint *> _neighbours; 

  list<double> data; //can be used to store extra-data attached to the point  

 private:
  Pt _coord;
  const int _no;
  float _value;
  
};

const bool operator ==(const Mpoint &p2, const Mpoint &p1);
const Vec operator -(const Mpoint&p1, const Mpoint &p2);
const bool operator <(const Mpoint &p1, const Mpoint &p2); //checks if p1 and p2 are adjacent
const Vec operator -(const Pt&p1, const Mpoint &p2);
const Vec operator -(const Mpoint&p1, const Pt &p2);

}

#endif
