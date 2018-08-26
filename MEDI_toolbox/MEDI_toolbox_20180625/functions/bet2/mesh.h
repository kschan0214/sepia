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

#ifndef _mesh
#define _mesh

#include <list>
#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include<algorithm>

#include "point.h"
#include "newimage/newimageall.h"

using namespace NEWIMAGE;

using namespace std;

namespace mesh{

class Mpoint;
class Triangle;

class Mesh {
 public:
  vector<Mpoint *> _points;
  list<Triangle *> _triangles;

  Mesh();
  ~Mesh();
  Mesh(const Mesh&m);
  Mesh operator=(const Mesh&m);

  void display() const;
  void clear();                 //clear the mesh and delete its components
  const int nvertices() const;
  Mpoint * get_point(int n){return _points[n];};
  
  double distance(const Pt& p) const; //signed distance of the point to the mesh
  void reorientate();     //puts the triangles in a coherent orientation
  void addvertex(Triangle *const t,const Pt p);
  void retessellate(); //global retesselation
  void update();       //puts _update_coord into _coords for each point  
  void translation(const double x,const double y,const double z);
  void translation(const Vec v);
  void rotation(const double r11, const double r12, const double r13,const double r21, const double r22, const double r23,const double r31, const double r32, const double r33, const double x, const double y, const double z);

  void rescale(const double t, const double x=0,const double y=0,const double z=0);
  void rescale(const double t , const Pt p);
  int load(string s="manual_input"); //  returns: -1 if load fails, 0 if load is cancelled
  //  1 if load succeeds and file is a .off file, 2 if load succeeds and file is a freesurfer file
  void load_off(string s="manual_input");
  void load_vtk_ASCII(string s="manual_input");
  void load_fs(string s="manual_input");
  void load_fs_label(string s="manual_input",const int& value=1);
  void save(string s="manual_input",int type=1) const;
  void save_fs_label(string s, bool RAS=false) const;//save an fs label of all points whos value greater than 0
                                                                        //If RAS is true, then mesh points are in standard-space coords
                                                                        // and these are saved as label coords
  const double self_intersection(const Mesh& original) const;
  const bool real_self_intersection();
  void stream_mesh(ostream& flot, int type=1) const; //type=1 -> .off style stream. type=2 -> freesurfer style stream
};

void make_mesh_from_tetra(int, Mesh&);
void make_mesh_from_icosa(int, Mesh&);
void make_mesh_from_octa(int, Mesh&);
ostream& operator <<(ostream& flot,const Mesh & m);

}

#endif


