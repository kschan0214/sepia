/*  Templated iterators for the image storage class (which uses lazymanager)

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 2000 University of Oxford  */

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

#if !defined(__positerators_h)
#define __positerators_h

#include <cstdlib>
#include "lazy.h"
#include <iostream>

using namespace LAZY;

namespace NEWIMAGE {

 //---------------------------------------------------------------------//

  // Mutable, Forward Iterator

  template <class T>
  class poslazyiterator {
  private:
    T* iter;
    lazymanager* lazyptr;
    int x;
    int y;
    int z;
    int lx0;
    int ly0;
    int lz0;
    int lx1;
    int ly1;
    int lz1;
    int RowAdjust;
    int SliceAdjust;

    void calc_offsets(int rowoff, int sliceoff) {
      RowAdjust = rowoff + (lx0 - lx1);
      SliceAdjust = sliceoff + rowoff*(ly0 - ly1) + (lx0 - lx1); 
    }

  public:
    poslazyiterator() : lazyptr(0) { }
    poslazyiterator(const poslazyiterator<T>& source)
      { this->operator=(source); }
    poslazyiterator(T* sourceptr, lazymanager* lazyp, 
		    int xinit, int yinit, int zinit,
		    int x0, int y0, int z0, int x1, int y1, int z1,
		    int rowoff, int sliceoff) :
      iter(sourceptr), lazyptr(lazyp), x(xinit), y(yinit), z(zinit),
      lx0(x0), ly0(y0), lz0(z0), lx1(x1), ly1(y1), lz1(z1)
      { calc_offsets(rowoff,sliceoff); }
    ~poslazyiterator() { }  // do nothing

    inline const poslazyiterator<T> operator++(int) 
      { poslazyiterator<T> tmp=*this; ++(*this); return tmp; }
    inline const poslazyiterator<T>& operator++() // prefix
      { x++; if (x>lx1) { x=lx0; y++; 
               if (y>ly1) { y=ly0; z++; if (z>lz1) { ++iter; } // end condition
                   else { iter+=SliceAdjust; } } 
               else { iter+=RowAdjust; } } 
             else { ++iter;}   return *this; }

    inline bool operator==(const poslazyiterator<T>& it) const
       { return iter == it.iter; }
    inline bool operator!=(const poslazyiterator<T>& it) const 
       { return iter != it.iter; }

    inline void getposition(int &rx, int &ry, int &rz) const
      { rx= x; ry = y; rz = z; }
    inline const int& getx() const { return x; }
    inline const int& gety() const { return y; }
    inline const int& getz() const { return z; }

    inline const poslazyiterator<T>& 
      operator=(const poslazyiterator<T>& source) 
      { iter = source.iter; lazyptr = source.lazyptr;
        lx0=source.lx0; ly0=source.ly0; lz0=source.lz0; 
	lx1=source.lx1; ly1=source.ly1; lz1=source.lz1; 
	x=source.x; y=source.y; z=source.z; 
	RowAdjust = source.RowAdjust;  SliceAdjust = source.SliceAdjust;
	return *this; }

    inline T& operator*() const 
      { lazyptr->set_whole_cache_validity(false); return *iter;}
  };


  //---------------------------------------------------------------------//

  // Constant, Forward Iterator

  template <class T>
  class posconstiterator {
  private:
    T* iter;
    int x;
    int y;
    int z;
    int lx0;
    int ly0;
    int lz0;
    int lx1;
    int ly1;
    int lz1;
    int RowAdjust;
    int SliceAdjust;

    void calc_offsets(int rowoff, int sliceoff) {
      RowAdjust = rowoff + (lx0 - lx1);
      SliceAdjust = sliceoff + rowoff*(ly0 - ly1) + (lx0 - lx1); 
    }

  public:
    posconstiterator() { }
    posconstiterator(const posconstiterator<T>& source)
      { this->operator=(source); }
    posconstiterator(T* sourceptr,
		    int xinit, int yinit, int zinit,
		    int x0, int y0, int z0, int x1, int y1, int z1,
		    int rowoff, int sliceoff) :
      iter(sourceptr), x(xinit), y(yinit), z(zinit),
      lx0(x0), ly0(y0), lz0(z0), lx1(x1), ly1(y1), lz1(z1)
      { calc_offsets(rowoff,sliceoff); }
    ~posconstiterator() { }  // do nothing

    inline const posconstiterator<T> operator++(int) 
      { posconstiterator<T> tmp=*this; ++(*this); return tmp; }
    inline const posconstiterator<T>& operator++() // prefix
      { x++; if (x>lx1) { x=lx0; y++; 
               if (y>ly1) { y=ly0; z++; if (z>lz1) { ++iter; } // end condition
                   else { iter+=SliceAdjust; } } 
               else { iter+=RowAdjust; } } 
             else { ++iter;}   return *this; }

    inline bool operator==(const posconstiterator<T>& it) const
       { return iter == it.iter; }
    inline bool operator!=(const posconstiterator<T>& it) const 
       { return iter != it.iter; }

    inline void getposition(int &rx, int &ry, int &rz) const
      { rx= x; ry = y; rz = z; }
    inline const int& getx() const { return x; }
    inline const int& gety() const { return y; }
    inline const int& getz() const { return z; }

    inline const posconstiterator<T>& 
      operator=(const posconstiterator<T>& source) 
      { iter = source.iter; 
        lx0=source.lx0; ly0=source.ly0; lz0=source.lz0; 
	lx1=source.lx1; ly1=source.ly1; lz1=source.lz1; 
	x=source.x; y=source.y; z=source.z; 
	RowAdjust = source.RowAdjust;  SliceAdjust = source.SliceAdjust;
	return *this; }

    inline const T& operator*() const { return *iter;}
  };



  //---------------------------------------------------------------------//

}  // end namespace

#endif

