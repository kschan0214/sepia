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

#include "mesh.h"
#include "triangle.h"
#include "mpoint.h"
#include "pt_special.h"


namespace mesh {

Mesh::Mesh(){
}

Mesh::~Mesh(){
  for (list<Triangle*>::iterator i=_triangles.begin(); i!=_triangles.end(); i++)
    delete (*i);
  for (vector<Mpoint*>::iterator i=_points.begin(); i!=_points.end(); i++)
    delete(*i);
}

Mesh::Mesh(const Mesh&m)
{
  _points.clear();
  _triangles.clear();
  for (vector<Mpoint*>::const_iterator p= m._points.begin(); p!=m._points.end(); p++)
    {
      Mpoint * pt = new Mpoint((*p)->get_coord(), (*p)->get_no());
      _points.push_back(pt);
    }
  for (list<Triangle*>::const_iterator t= m._triangles.begin(); t!=m._triangles.end(); t++)
    {
      int v0 = (*t)->get_vertice(0)->get_no(), v1 = (*t)->get_vertice(1)->get_no(), v2 = (*t)->get_vertice(2)->get_no();
      Triangle * tr = new Triangle(get_point(v0), get_point(v1), get_point(v2));
      _triangles.push_back(tr);
    }
}

 

Mesh Mesh::operator=(const Mesh&m)
{
  if (this == &m) return *this;
  for (list<Triangle*>::iterator i=_triangles.begin(); i!=_triangles.end(); i++)
    delete (*i);
  for (vector<Mpoint*>::iterator i=_points.begin(); i!=_points.end(); i++)
    delete(*i);
  _points.clear();
  _triangles.clear();
  for (vector<Mpoint*>::const_iterator p= m._points.begin(); p!=m._points.end(); p++)
    {
      Mpoint * pt = new Mpoint((*p)->get_coord(), (*p)->get_no());
      (*pt).data=(*p)->data; //optional line to copy data
      _points.push_back(pt);
    }
  for (list<Triangle*>::const_iterator t= m._triangles.begin(); t!=m._triangles.end(); t++)
    {
      int v0 = (*t)->get_vertice(0)->get_no(), v1 = (*t)->get_vertice(1)->get_no(), v2 = (*t)->get_vertice(2)->get_no();
      Triangle * tr = new Triangle(get_point(v0), get_point(v1), get_point(v2));
      _triangles.push_back(tr);
    }
  
  return *this;
}


const int Mesh::nvertices() const {
  return _points.size();
}



void Mesh::display() const{
  cout<<*this<<endl;
}

void Mesh::clear() {
  for (list<Triangle*>::iterator i=_triangles.begin(); i!=_triangles.end(); i++)
    delete (*i);
  _triangles.clear();
  for (vector<Mpoint*>::iterator i=_points.begin(); i!=_points.end(); i++)
    delete(*i);
  _points.clear();
}



void make_mesh_from_octa(int n, Mesh& m)
{
  m.clear();
  
  Mpoint *XPLUS = new Mpoint( 1,  0,  0, 0);
  Mpoint *XMIN  = new Mpoint(-1,  0,  0, 1);
  Mpoint *YPLUS = new Mpoint( 0,  1,  0, 2);
  Mpoint *YMIN  = new Mpoint( 0, -1,  0, 3);
  Mpoint *ZPLUS = new Mpoint( 0,  0,  1, 4);
  Mpoint *ZMIN  = new Mpoint( 0,  0, -1, 5);

  Triangle * t0= new Triangle( XPLUS, ZPLUS, YPLUS );
  Triangle * t1= new Triangle( YPLUS, ZPLUS, XMIN  );
  Triangle * t2= new Triangle( XMIN , ZPLUS, YMIN  );
  Triangle * t3= new Triangle( YMIN , ZPLUS, XPLUS );
  Triangle * t4= new Triangle( XPLUS, YPLUS, ZMIN  );
  Triangle * t5= new Triangle( YPLUS, XMIN , ZMIN  );
  Triangle * t6= new Triangle( XMIN , YMIN , ZMIN  );
  Triangle * t7= new Triangle( YMIN , XPLUS, ZMIN  );

  m._points.push_back(XPLUS);
  m._points.push_back(XMIN);
  m._points.push_back(YPLUS);
  m._points.push_back(YMIN);
  m._points.push_back(ZPLUS);
  m._points.push_back(ZMIN);
  
  m._triangles.push_back(t0);
  m._triangles.push_back(t1);
  m._triangles.push_back(t2);
  m._triangles.push_back(t3);
  m._triangles.push_back(t4);
  m._triangles.push_back(t5);
  m._triangles.push_back(t6);
  m._triangles.push_back(t7);
   
  for (int i = 1; i<n ; i++)
    {
      //re-tesselates
      m.retessellate();
      
      for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++)
	{
	  double norm = sqrt((*i)->get_coord().X * (*i)->get_coord().X + (*i)->get_coord().Y * (*i)->get_coord().Y + (*i)->get_coord().Z * (*i)->get_coord().Z);
	  (*i)->_update_coord=(*i)->get_coord();
	  (*i)->_update_coord*=(1/norm);
	}
      m.update();
    }
}


void make_mesh_from_icosa(int n, Mesh& m)
{
  m.clear();
  
  const double tau=0.8506508084;	
  const double one=0.5257311121;
  Mpoint *ZA=new Mpoint( tau , one,  0,    0);
  Mpoint *ZB=new Mpoint( -tau, one,  0,    1);
  Mpoint *ZC=new Mpoint( -tau, -one, 0,    2);
  Mpoint *ZD=new Mpoint( tau,  -one, 0,    3);
  Mpoint *YA=new Mpoint( one,  0 ,   tau,  4);
  Mpoint *YB=new Mpoint( one,  0,    -tau, 5);
  Mpoint *YC=new Mpoint( -one, 0,    -tau, 6);
  Mpoint *YD=new Mpoint( -one, 0,    tau,  7);
  Mpoint *XA=new Mpoint( 0 ,   tau,  one,  8);
  Mpoint *XB=new Mpoint( 0,    -tau, one,  9);
  Mpoint *XC=new Mpoint( 0,    -tau, -one, 10);
  Mpoint *XD=new Mpoint( 0,    tau,  -one, 11);

  Triangle * t0= new Triangle(YA, XA, YD);
  Triangle * t1= new Triangle(YA, YD, XB);
  Triangle * t2= new Triangle(YB, YC, XD);
  Triangle * t3= new Triangle(YB, XC, YC);
  Triangle * t4= new Triangle(ZA, YA, ZD);
  Triangle * t5= new Triangle(ZA, ZD, YB);
  Triangle * t6= new Triangle(ZC, YD, ZB);
  Triangle * t7= new Triangle(ZC, ZB, YC);
  Triangle * t8= new Triangle(XA, ZA, XD);
  Triangle * t9= new Triangle(XA, XD, ZB);
  Triangle * t10= new Triangle(XB, XC, ZD);
  Triangle * t11= new Triangle(XB, ZC, XC);
  Triangle * t12= new Triangle(XA, YA, ZA);
  Triangle * t13= new Triangle(XD, ZA, YB);
  Triangle * t14= new Triangle(YA, XB, ZD);
  Triangle * t15= new Triangle(YB, ZD, XC);
  Triangle * t16= new Triangle(YD, XA, ZB);
  Triangle * t17= new Triangle(YC, ZB, XD);
  Triangle * t18= new Triangle(YD, ZC, XB);
  Triangle * t19= new Triangle(YC, XC, ZC);

  m._points.push_back(ZA);
  m._points.push_back(ZB);
  m._points.push_back(ZC);
  m._points.push_back(ZD);
  m._points.push_back(YA);
  m._points.push_back(YB);
  m._points.push_back(YC);
  m._points.push_back(YD);
  m._points.push_back(XA);
  m._points.push_back(XB);
  m._points.push_back(XC);
  m._points.push_back(XD);

  m._triangles.push_back(t0);
  m._triangles.push_back(t1);
  m._triangles.push_back(t2);
  m._triangles.push_back(t3);
  m._triangles.push_back(t4);
  m._triangles.push_back(t5);
  m._triangles.push_back(t6);
  m._triangles.push_back(t7);
  m._triangles.push_back(t8);
  m._triangles.push_back(t9);
  m._triangles.push_back(t10);
  m._triangles.push_back(t11);
  m._triangles.push_back(t12);
  m._triangles.push_back(t13);
  m._triangles.push_back(t14);
  m._triangles.push_back(t15);
  m._triangles.push_back(t16);
  m._triangles.push_back(t17);
  m._triangles.push_back(t18);
  m._triangles.push_back(t19);

  for (list<Triangle *>::iterator i=m._triangles.begin(); i!=m._triangles.end(); i++)
    (*i)->swap();

  for (int io = 1; io<n ; io++)
    {
      //re-tesselates
      m.retessellate();
      
      for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++)
	{
	  double norm = sqrt((*i)->get_coord().X * (*i)->get_coord().X + (*i)->get_coord().Y * (*i)->get_coord().Y + (*i)->get_coord().Z * (*i)->get_coord().Z);
	  (*i)->_update_coord=(*i)->get_coord();
	  (*i)->_update_coord*=(1/norm);
	}
      m.update();
    }
}


void make_mesh_from_tetra(int n, Mesh& m)
{
  const double sqrt_3=0.5773502692;
  const Pt PPP(  sqrt_3,  sqrt_3,  sqrt_3);
  const Pt MMP(  -sqrt_3,  -sqrt_3,  sqrt_3);
  const Pt MPM(  -sqrt_3,  sqrt_3,  -sqrt_3);
  const Pt PMM(  sqrt_3,  -sqrt_3,  -sqrt_3);
  
  m.clear();
  
  Mpoint *p0=new Mpoint(PPP, 0);
  Mpoint *p1=new Mpoint(MMP, 1);
  Mpoint *p2=new Mpoint(MPM, 2);
  Mpoint *p3=new Mpoint(PMM, 3);
  Triangle * t0= new Triangle(p0, p1, p2);
  Triangle * t1= new Triangle(p0, p1, p3);
  Triangle * t2= new Triangle(p0, p3, p2);
  Triangle * t3= new Triangle(p3, p1, p2);
  m._points.push_back(p0);
  m._points.push_back(p1);
  m._points.push_back(p2);
  m._points.push_back(p3);
  m._triangles.push_back(t0);
  m._triangles.push_back(t1);
  m._triangles.push_back(t2);
  m._triangles.push_back(t3);
  
  
  for (int i = 1; i<n ; i++)
    {
      //re-tesselates
      m.retessellate();
      
      for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++)
	{
	  double norm = sqrt((*i)->get_coord().X * (*i)->get_coord().X + (*i)->get_coord().Y * (*i)->get_coord().Y + (*i)->get_coord().Z * (*i)->get_coord().Z);
	  (*i)->_update_coord=(*i)->get_coord();
	  (*i)->_update_coord*=(1/norm);
	}
      m.update();
    }
}

int Mesh::load(string s) {
  //loads a mesh - 
  //  returns: 
  // -1 if load fails, 
  // 0 if load is cancelled
  //  1 if load succeeds and file is a .off file, 
  //  2 if load succeeds and file is a freesurfer file
  // 3 if load succeed and file is a .vtk file
  clear();
  if (s == "manual_input")
    {
      // reads a line in standard input
      cout << "loading mesh : enter file name / c to cancel: ";
      s="";
      while (s.empty())
	{
	  string input;
	  getline(cin, input);
	  s = input;
	}
    }  
  //find out if it is an off file
  bool is_off=true;
  bool is_vtk=false;
  int ret=0;
  if (s!="c")
    {
      ifstream f(s.c_str());
      if (f.is_open())
	{	
	  //reading the header
	  string header;
	  getline(f, header);
	  //cout<<header;
	  string::size_type pos = header.find("OFF");
	  if (pos == string::npos) {
	    is_off=false;
	    pos=header.find("# vtk DataFile Version 3.0");
	    is_vtk=true;
	    if (pos == string::npos) {
	      is_vtk=false;
	    } 	  
	    
	  }
	  
	  f.close();
	  if(is_off){
	    //cout<<"Reading OFF file";
	    load_off(s);
	    ret=1;
	  }
	  else if(is_vtk){
	    //cout<<"Reading vtk file";
	    load_vtk_ASCII(s);
	    ret=3;
	  }
	  
	  else{
	    cout<<"Read other";
	    load_fs(s);
	    ret=2;
	  }
	}
      else {cout<<"error opening file"<<endl; ret=-1; exit(-1);}
    }

  else{ret=0; cout<<"cancelled"<<endl; }
  return ret;
}

void Mesh::load_off(string s) {//loads a .off format mesh
  clear();
  if (s == "manual_input")
    {
      // reads a line in standard input
      cout << "loading mesh : enter file name / c to cancel: ";
      s="";
      while (s.empty())
	{
	  string input;
	  getline(cin, input);
	  s = input;
	}
    }  

  // reads the data
  if (s!="c")
    {
      ifstream f(s.c_str());
      if (f.is_open())
	{	
	  //reading the header
	  string header;
	  getline(f, header);
	  string::size_type pos = header.find("OFF");
	  if (pos == string::npos) {cerr<<"error in the header"<<endl;exit(-1);}; 
	  pos = header.find("n");
	  if (pos != string::npos) {int N; f>>N; if (N!=3) {cerr<<"this program only handles triangles meshes"<<endl; exit(-1);}};

	  //reading the size of the mesh
	  int NVertices, NFaces, NEdges = 0;
	  f>>NVertices>>NFaces>>NEdges;

	  //reading the points
	  for (int i=0; i<NVertices; i++)
	    {
	      double x, y, z;
	      f>>x>>y>>z;
	      Mpoint * m = new Mpoint(x, y, z, i);
	      _points.push_back(m);
	    }
	  
	  //reading the triangles
	  for (int i=0; i<NFaces; i++)
	    {
	      int p0, p1, p2;
	      int j;
	      f>>j>>p0>>p1>>p2;
	      Triangle * t = new Triangle(get_point(p0), get_point(p1), get_point(p2));
	      _triangles.push_back(t);
	    }
	  f.close();
	}
      else {cout<<"error opening file"<<endl; exit(-1);}
    }
  else cout<<"cancelled"<<endl; 
}

void Mesh::load_vtk_ASCII(string s) {//loads a .vtk format mesh
  clear();
  if (s == "manual_input")
    {
      // reads a line in standard input
      cout << "loading mesh : enter file name / c to cancel: ";
      s="";
      while (s.empty())
	{
	  string input;
	  getline(cin, input);
	  s = input;
	}
    }  

  // reads the data
  if (s!="c")
    {
      ifstream f(s.c_str());
      if (f.is_open())
	{	
	  //reading the header
	  string header;
	  getline(f, header);
	  string::size_type pos = header.find("vtk DataFile Version 3.0");
	  if (pos == string::npos) {cerr<<"error in the header"<<endl;exit(-1);}
	  getline(f,header);
	  getline(f,header);
	  getline(f,header);
	  //getline(f,header);
	  //  cout<<header<<endl;
	  //reading the size of the mesh
	  int NVertices, NFaces;
	  f>>header>>NVertices>>header;
	  // cout<<"load vtkmesh npts:"<<NVertices;
	  
	  //reading the points
	  for (int i=0; i<NVertices; i++)
	    {
	      double x, y, z;
	      f>>x>>y>>z;
	      Mpoint * m = new Mpoint(x, y, z, i);
	      _points.push_back(m);
	    }
	  f>>header>>NFaces>>header;
	  //reading the triangles
	  for (int i=0; i<NFaces; i++)
	    {
	      int p0, p1, p2;
	      int j;
	      f>>j>>p0>>p1>>p2;
	      Triangle * t = new Triangle(get_point(p0), get_point(p1), get_point(p2));
	      _triangles.push_back(t);
	    }
	  f.close();
	}
      else {cout<<"error opening file"<<endl; exit(-1);}
    }
  else cout<<"cancelled"<<endl; 
}

void Mesh::load_fs(string s) { //load a freesurfer ascii mesh
  clear();

  if (s == "manual_input")
    {
      // reads a line in standard input
      cout << "loading mesh : enter file name / c to cancel: ";
      s="";
      while (s.empty())
	{
	  string input;
	  getline(cin, input);
	  s = input;
	}
    }  
  // reads the data
  if (s!="c")
    {
      ifstream f(s.c_str());
      if (f.is_open())
	{	
	  //reading the header
	  string header;
	  getline(f, header);
	  // string::size_type pos = header.find("OFF");
	  // if (pos == string::npos) {cerr<<"error in the header"<<endl;exit(-1);}; 
	  // pos = header.find("n");
	  // if (pos != string::npos) {int N; f>>N; if (N!=3) {cerr<<"this program only handles triangles meshes"<<endl; exit(-1);}};

	  //reading the size of the mesh
	  int NVertices, NFaces;
	  f>>NVertices>>NFaces;

	  //reading the points
	  for (int i=0; i<NVertices; i++)
	    {
	     double x, y, z;
	      float val;
	      f>>x>>y>>z>>val;
	      Mpoint * m = new Mpoint(x, y, z, i, val);
	      _points.push_back(m);
	    }
	  
	  //reading the triangles
	  for (int i=0; i<NFaces; i++)
	    {
	      int p0, p1, p2;
	      float val;
	      //	      int j;
	      f>>p0>>p1>>p2>>val;
	      Triangle * t = new Triangle(get_point(p0), get_point(p1), get_point(p2),val);
	      _triangles.push_back(t);
	    }
	  f.close();
	}
      else {cout<<"error opening file"<<endl; exit(-1);}
    }
  else cout<<"cancelled"<<endl; 
}


void Mesh::load_fs_label(string s,const int& value){
  if (s == "manual_input")
    {
      // reads a line in standard input
      cout << "loading label : enter file name / c to cancel: ";
      s="";
      while (s.empty())
	{
	  string input;
	  getline(cin, input);
	  s = input;
	}
    }  
   // reads the data
  if (s!="c")
    {
      ifstream f(s.c_str());
      if (f.is_open())
	{	
	  //reading the header
	  string header;
	  getline(f, header);
	  // string::size_type pos = header.find("OFF");
	  // if (pos == string::npos) {cerr<<"error in the header"<<endl;exit(-1);}; 
	  // pos = header.find("n");
	  // if (pos != string::npos) {int N; f>>N; if (N!=3) {cerr<<"this program only handles triangles meshes"<<endl; exit(-1);}};

	  //reading the size of the mesh
	  int NVertices;
	  f>>NVertices;

	  //reading the points
	  for (int i=0; i<NVertices; i++)
	    {
	      int num;
	      double x, y, z;
	      float tmp;
	      f>>num>>x>>y>>z>>tmp; //NB - can't work out when Freesurfer sets the 5th value in label
	      _points[num]->set_value(value);
	    }
	  
	  f.close();
	}
      else {cout<<"error opening file"<<endl; exit(-1);}
    }
  else cout<<"cancelled"<<endl; 
}


void Mesh::save(string s,int type) const {//type is 1 for an off file, 2 for a FS file
  if (strcmp(s.c_str(), "c"))
    {
      // writes the file
      ofstream f(s.c_str());
      if (f.is_open())
	{
	  stream_mesh(f,type);
	  f.close();
	}
      else cerr<<"error opening file "<<s<<endl;
    }
  else cerr<<"cancelled"<<endl;
  
}

void Mesh::save_fs_label(string s, bool RAS) const { //save an fs label of all points whos value greater than 0
                                                                      //If RAS is true, then mesh points are in standard-space coords
                                                                      // and these are saved as label coords
  ofstream f(s.c_str());
  stringstream flot;
  if (f.is_open())
    {
      int ptcount=0,labcount=0;
      for(vector<Mpoint *>::const_iterator i=_points.begin();i!=_points.end();i++){
	if((*i)->get_value()>0){
	  if(!RAS)
	    flot<<ptcount<<" "<<0<<" "<<0<<" "<<0<<" "<<(*i)->get_value()<<endl; //Don't know RAS coordinate of point
	  else
	    flot<<ptcount<<" "<<(*i)->get_coord().X<<" "<<(*i)->get_coord().Y<<" "<<(*i)->get_coord().Y<<" "<<(*i)->get_value()<<endl; 
	  
	  labcount++;
	}
	
	ptcount++;
	
	
      }
      f<<"#!ascii label , from subject"<<endl;
      f<<labcount<<endl<<flot.str();      
      f.close();
    }
  else cerr<<"error opening file "<<s<<endl;
}

void Mesh::addvertex(Triangle * const t,const Pt p)
{
  Mpoint *_p;
  _p = new Mpoint(p, nvertices());
  Triangle *_t1, *_t2, *_t3;
  _t1 = new Triangle(_p, t->get_vertice(1), t->get_vertice(0));
  _t2 = new Triangle(_p, t->get_vertice(0), t->get_vertice(2));
  _t3 = new Triangle(_p, t->get_vertice(2), t->get_vertice(1));

  _triangles.remove(t);
  delete(t);

  _points.push_back(_p);
  _triangles.push_back(_t1);
  _triangles.push_back(_t2);
  _triangles.push_back(_t3);
}

//gives a correct orientation to the triangles
void Mesh::reorientate() {
  int count = 0;
  list<Triangle *> prov = _triangles;
  while (!prov.empty())
    {
      count++;
      Triangle * current = prov.front();
      prov.remove(current);
      current->oriented = true;
      if(!(prov.empty()))
	for (int c=0; c<3; c++)
	  for (list<Triangle*>::iterator i = current->get_vertice(c)->_triangles.begin(); i!= current->get_vertice(c)->_triangles.end();i++ )
	    {
	      int m;
	      {m = *(*i)<current;}
	      switch (m) {
	      case 0: 
		break;
	      case 1: 
		if (!(*i)->oriented) prov.push_front(*i);
		break;
	      case 2:
		if (!(*i)->oriented) {(*i)->swap();
		prov.push_front(*i);}
		break;
	      }
	    }
    }
  for (list<Triangle*>::iterator i = _triangles.begin(); i!= _triangles.end();i++ )
    (*i)->oriented = false;
}

void Mesh::update()
{
  for (vector<Mpoint *>::iterator i = _points.begin(); i!=_points.end(); i++)
    (*i)->update();
}

void Mesh::translation(const double x,const double y,const double z)
{
  for (vector<Mpoint *>::iterator i = _points.begin(); i!=_points.end(); i++)
    (*i)->translation(x, y, z);
}

void Mesh::translation(const Vec v)
{
  double x = v.X, y = v.Y, z = v.Z;
  for (vector<Mpoint *>::iterator i = _points.begin(); i!=_points.end(); i++)
    (*i)->translation(x, y, z);
}

void Mesh::rotation(const double r11, const double r12, const double r13,const double r21, const double r22, const double r23,const double r31, const double r32, const double r33, const double x, const double y, const double z)
{
  for (vector<Mpoint *>::iterator i = _points.begin(); i!=_points.end(); i++)
    (*i)->rotation(r11,r12,r13,r21,r22,r23,r31,r32,r33, x, y, z);
}

void Mesh::rescale(const double t, const double x, const double y, const double z)
{
  for (vector<Mpoint *>::iterator i = _points.begin(); i!=_points.end(); i++)
    (*i)->rescale(t, x, y, z);
}

void Mesh::rescale(const double t , const Pt p)
{
  for (vector<Mpoint *>::iterator i = _points.begin(); i!=_points.end(); i++)
    (*i)->rescale(t, p.X, p.Y, p.Z);
}

void Mesh::retessellate() {

  vector<Mpoint *> added_points;
  list<Triangle*> tr = _triangles;

  added_points.clear();

  for (list<Triangle*>::iterator t=tr.begin(); t!=tr.end(); t++)
    {
      Mpoint* v0=(*t)->get_vertice(0);
      Mpoint* v1=(*t)->get_vertice(1);
      Mpoint* v2=(*t)->get_vertice(2);
      Pt pt0((v1->get_coord().X + v2->get_coord().X)/2, 
	     (v1->get_coord().Y + v2->get_coord().Y)/2,
	     (v1->get_coord().Z + v2->get_coord().Z)/2);
      Pt pt2((v0->get_coord().X + v1->get_coord().X)/2, 
	     (v0->get_coord().Y + v1->get_coord().Y)/2,
	     (v0->get_coord().Z + v1->get_coord().Z)/2);
      Pt pt1((v0->get_coord().X + v2->get_coord().X)/2, 
	     (v0->get_coord().Y + v2->get_coord().Y)/2,
	     (v0->get_coord().Z + v2->get_coord().Z)/2);

      Mpoint* p1 = NULL;
      Mpoint* p2 = NULL;
      Mpoint* p0 = NULL;
      
      int count=0;

      bool b0=true, b1=true, b2=true;

      for (vector<Mpoint*>::const_iterator mi=added_points.begin(); mi!=added_points.end(); mi++)
	{
	  Pt current = (*mi)->get_coord();
	  if (pt0==current) {b0=false; p0=*mi;}
	  if (pt1==current) {b1=false; p1=*mi;}
	  if (pt2==current) {b2=false; p2=*mi;}
	}


      if (b0) {p0=new Mpoint(pt0, nvertices() + count); count++;};
      if (b1) {p1=new Mpoint(pt1, nvertices() + count); count++;};
      if (b2) {p2=new Mpoint(pt2, nvertices() + count); count++;}; 
      
      Triangle* t0 = new Triangle(p2, p0, p1);
      Triangle* t1 = new Triangle(p1, v0, p2);
      Triangle* t2 = new Triangle(p0, v2, p1);
      Triangle* t3 = new Triangle(p2, v1, p0);
      
      _triangles.push_back(t0);
      _triangles.push_back(t1);
      _triangles.push_back(t2);
      _triangles.push_back(t3);
      
      if (b0) {_points.push_back(p0);added_points.push_back(p0);}
      if (b1) {_points.push_back(p1);added_points.push_back(p1);}
      if (b2) {_points.push_back(p2);added_points.push_back(p2);}
      
      v0->_neighbours.remove(v1);
      v0->_neighbours.remove(v2);
      v1->_neighbours.remove(v0);
      v1->_neighbours.remove(v2);
      v2->_neighbours.remove(v1);
      v2->_neighbours.remove(v0);

      vector<Mpoint*>::iterator ite2=added_points.begin();
      ite2++;
    }

  for (list<Triangle*>::iterator t=tr.begin(); t!=tr.end(); t++)
    {
      _triangles.remove(*t);
      delete (*t);
    }
}


double Mesh::distance(const Pt& p) const
{
  bool triangle = false;
  Mpoint * nearest_point = NULL;
  Triangle * nearest_triangle = NULL;
  double min = 50000; // ...
  
  for (vector<Mpoint *>::const_iterator i = _points.begin(); i != _points.end(); i++)
    if (((**i) - p).norm() < min)
      {
	min = ((**i) - p).norm();
	nearest_point = (*i);
      }
  
  for (list<Triangle *>::const_iterator i = _triangles.begin(); i!=_triangles.end(); i++)
    {
      Pt pp; //projection of p on the triangle
      double d = 50000;
      double a, b, c, xm, ym, zm;
      Vec n = (*i)->normal(); n.normalize();
      a=n.X; b=n.Y; c=n.Z;
      Pt M = (*i)->get_vertice(0)->get_coord(); xm=M.X; ym=M.Y; zm=M.Z;
      pp.Z =p.Z + c* (a* (xm - p.X) + b* (ym - p.Y) + c* (zm - p.Z));
      pp.X =p.X + a* (a* (xm - p.X) + b* (ym - p.Y) + c* (zm - p.Z));
      pp.Y =p.Y + b* (a* (xm - p.X) + b* (ym - p.Y) + c* (zm - p.Z));
      
      Vec n0, n1, n2;
      n0 = (*(*i)->get_vertice(2) - *(*i)->get_vertice(1))*(*(*i)->get_vertice(2) - pp);
      n1 = (*(*i)->get_vertice(0) - *(*i)->get_vertice(2))*(*(*i)->get_vertice(0) - pp);
      n2 = (*(*i)->get_vertice(1) - *(*i)->get_vertice(0))*(*(*i)->get_vertice(1) - pp);
      if ((n0|n1) >= 0 && (n0|n2) >= 0) d = (pp - p).norm();
      
      if (d<min)
	{
	  min = d;
	  nearest_triangle = (*i);
	  triangle = true;
	}
      
    }
 
  if (triangle) {if ((nearest_triangle->normal()|(nearest_triangle->centroid() - p))>0) min=-min;}
  else {if ((nearest_point->local_normal()|(*nearest_point - p))>0) min=-min;}
  
  
  return min;
}


const double Mesh::self_intersection(const Mesh& original) const
{
  double intersection = 0;
  if (original._points.size() != _points.size()) {cerr<<"error, parameter for self_intersection should be the original mesh"<<endl; return -1;};
  
  vector<Mpoint*>::const_iterator io = original._points.begin();
  
  double ml = 0, mlo = 0;
  int counter = 0;
  for (vector<Mpoint*>::const_iterator i = _points.begin(); i!=_points.end(); i++, io++ )
    {
      counter++;
      ml += (*i)->medium_distance_of_neighbours();
      mlo += (*io)->medium_distance_of_neighbours();
    }
  ml/=counter;
  mlo/=counter;


  io = original._points.begin();
  vector<Mpoint*>::const_iterator jo = original._points.begin();
  
   for (vector<Mpoint*>::const_iterator i = _points.begin(); i!=_points.end(); i++, io++)
     {
       jo = original._points.begin();
    for (vector<Mpoint*>::const_iterator j = _points.begin(); j!=_points.end(); j++, jo++)
      if (!((*i)==(*j)) && !((**i)<(**j)))
	if (((*i)->get_coord().X -  (*j)->get_coord().X) * ((*i)->get_coord().X -  (*j)->get_coord().X) + ((*i)->get_coord().Y -  (*j)->get_coord().Y) * ((*i)->get_coord().Y -  (*j)->get_coord().Y) + ((*i)->get_coord().Z -  (*j)->get_coord().Z) * ((*i)->get_coord().Z -  (*j)->get_coord().Z)< ml * ml)
	  {
	    double dist = (((**i) - (**j)).norm())/ml;
	    double disto = (((**io) - (**jo)).norm())/mlo;
	    intersection += (dist - disto)*(dist - disto);
	  }
     }
   return intersection;
}

void Mesh::stream_mesh(ostream& flot,int type) const{ 
	flot.setf(ios::fixed);// 6 decimal places - matches freesurfer..
	flot.precision(6);
#ifdef PPC64
    int n=0;
#endif
	
	if(type==1){
		flot<<"OFF"<<endl;
		flot<<_points.size()<<" "<<_triangles.size()<<" "<<"0"<<endl;
		for (vector<Mpoint *>::const_iterator i =_points.begin(); i!=_points.end(); i++)  
		{ 
			//	flot.precision(6);
			flot<<(*i)->get_coord().X<<" "<<(*i)->get_coord().Y<<" "<<(*i)->get_coord().Z<<endl;
#ifdef PPC64
			if ((n++ % 20) == 0) flot.flush();
#endif
		}
		for ( list<Triangle*>::const_iterator i=_triangles.begin(); i!=_triangles.end(); i++) 
			flot<<"3 "<<(*i)->get_vertice(0)->get_no()<<" "<<(*i)->get_vertice(1)->get_no()<<" "<<(*i)->get_vertice(2)->get_no()<<" "<<endl;
#ifdef PPC64
        if ((n++ % 20) == 0) flot.flush();
#endif
	}
	else if(type==2){
		flot<<"#ascii FS mesh"<<endl;
		flot<<_points.size()<<" "<<_triangles.size()<<endl;
		for (vector<Mpoint *>::const_iterator i =_points.begin(); i!=_points.end(); i++)  
		{ 
			//	flot.precision(6);
			flot<<(*i)->get_coord().X<<" "<<(*i)->get_coord().Y<<" "<<(*i)->get_coord().Z<<" "<<(*i)->get_value()<<endl;
#ifdef PPC64
			if ((n++ % 20) == 0) flot.flush();
#endif
		}
		for ( list<Triangle*>::const_iterator i=_triangles.begin(); i!=_triangles.end(); i++) 
			flot<<(*i)->get_vertice(0)->get_no()<<" "<<(*i)->get_vertice(1)->get_no()<<" "<<(*i)->get_vertice(2)->get_no()<<" "<<(*i)->get_value()<<endl;
#ifdef PPC64
        if ((n++ % 20) == 0) flot.flush();
#endif
	} else if(type==3){
		flot<<"# vtk DataFile Version 3.0"<<endl<<"surface file"<<endl<<"ASCII"<<endl<<"DATASET POLYDATA"<<endl<<"POINTS ";
		flot<<_points.size()<<"  float"<<endl;
		
		for (vector<Mpoint *>::const_iterator i =_points.begin(); i!=_points.end(); i++)  
		{ 
			//	flot.precision(6);
			flot<<(*i)->get_coord().X<<" "<<(*i)->get_coord().Y<<" "<<(*i)->get_coord().Z<<endl;
#ifdef PPC64
			if ((n++ % 20) == 0) flot.flush();
#endif
		}
		flot<<"POLYGONS "<<_triangles.size()<<" "<<_triangles.size()*4<<endl;
		for ( list<Triangle*>::const_iterator i=_triangles.begin(); i!=_triangles.end(); i++) 
			flot<<"3 "<<(*i)->get_vertice(0)->get_no()<<" "<<(*i)->get_vertice(1)->get_no()<<" "<<(*i)->get_vertice(2)->get_no()<<" "<<endl;
#ifdef PPC64
        if ((n++ % 20) == 0) flot.flush();
#endif
	}
	
	else{cout<<"Invalid file Type"<<endl;exit(-1);}
}



ostream& operator <<(ostream& flot,const Mesh & m){
  m.stream_mesh(flot,1);
  return flot;
}

const bool Mesh::real_self_intersection()
{
  vector<Pt_special *> v;
  for (vector<Mpoint *>::const_iterator i = this->_points.begin(); i != this->_points.end(); i++)
    {
      Pt_special * p = new (Pt_special); 
      p->P = (*i);
      for (list<Triangle *>::const_iterator i2 = (*i)->_triangles.begin(); i2 != (*i)->_triangles.end(); i2++)
	{p->T.push_back(*i2);}
      v.push_back(p); 
    }
  


  for (list<Triangle *>::const_iterator i = this->_triangles.begin(); i != this->_triangles.end(); i++)
    {
      //initialize
      (*i)->data.clear();
      (*i)->data.push_back(0);

      //compute the center
      Pt centre(0, 0, 0);
      centre += (*i)->get_vertice(0)->get_coord();
      centre += (*i)->get_vertice(1)->get_coord();
      centre += (*i)->get_vertice(2)->get_coord();
      centre *= 1/3;
      (*i)->data.push_back(centre.X);
      (*i)->data.push_back(centre.Y);
      (*i)->data.push_back(centre.Z);
      
      //compute the max distance to the center
      double max = 0;
      for (int e = 0; e < 3; e++)
	{
	  double c = (centre.X - (*i)->get_vertice(e)->get_coord().X)*(centre.X - (*i)->get_vertice(e)->get_coord().X) + (centre.Y - (*i)->get_vertice(e)->get_coord().Y)*(centre.Y - (*i)->get_vertice(e)->get_coord().Y) + (centre.Z - (*i)->get_vertice(e)->get_coord().Z)*(centre.Z - (*i)->get_vertice(e)->get_coord().Z);
	  if (c > max) max = c;
	}
      (*i)->data.push_back(max);
      
    }
  sort(v.begin(), v.end(), compPt());
  
  bool result = false;
  
  //swipping
  list<Triangle *> l;
  l.clear();
  for (vector<Pt_special *>::iterator i = v.begin(); i != v.end(); i++)
    {
      if (result) break;
      for (list<Triangle *>::iterator tr = (*i)->T.begin(); tr != (*i)->T.end(); tr++)
	{
	  if (result) break;
	  if ((*tr)->data[0] == 0)
	    {
	      //main loop
	      for (list<Triangle *>::const_iterator tr2 = l.begin(); tr2 != l.end(); tr2++)
		{
		  bool res = false;
		  for (int i = 0; i < 3; i++)
		    for (int j = 0; j < 3; j++)
		      {
			if (!res) if ((*tr2)->get_vertice(i)->get_coord() == (*tr)->get_vertice(j)->get_coord()) res = true;
		      }
		  if (!res)
		    {
		      double d = ((*tr)->data[1] - (*tr2)->data[1]) * ((*tr)->data[1] - (*tr2)->data[1]) + ((*tr)->data[2] - (*tr2)->data[2]) * ((*tr)->data[2] - (*tr2)->data[2]) + ((*tr)->data[3] - (*tr2)->data[3]) * ((*tr)->data[3] - (*tr2)->data[3]);
		      if (d < (*tr)->data[4] + (*tr2)->data[4] )
			
			{
			  result = result | (*tr)->intersect(**tr2);
			}
		    }
		}
	      l.push_back(*tr);
	      (*tr)->data[0] = 1;
	    }
	  else 
	    if ((*tr)->data[0] == 1) (*tr)->data[0] = 2; 
	    else if ((*tr)->data[0] == 2) {l.remove(*tr);}
	}
    }

  return result;

}

}












