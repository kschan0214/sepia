/*  Templated iterators for any storage class that uses lazymanager

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

#if !defined(__lazyiterators_h)
#define __lazyiterators_h

#include <cstdlib>
#include "lazy.h"
#include <iostream>

using namespace LAZY;

namespace LAZY {

  // Mutable, Random Access Iterator

  template <class IT, class T>
  class rlazyiterator {
  private:
    IT iter;
    lazymanager* lazyptr;
  public:
    rlazyiterator() : lazyptr(0) { }
    rlazyiterator(const rlazyiterator<IT,T>& source) : 
      iter(source.iter) , lazyptr(source.lazyptr) { }
    rlazyiterator(const IT& sourceiter, lazymanager* lazyp) :
      iter(sourceiter) , lazyptr(lazyp) { }
    ~rlazyiterator() { }  // do nothing

    inline const rlazyiterator<IT,T> operator++(int) 
      { rlazyiterator<IT,T> tmp=*this; iter++; return tmp; }
    inline const rlazyiterator<IT,T>& operator++() // prefix
      { ++iter; return *this; }
    inline const rlazyiterator<IT,T> operator--(int) 
      { rlazyiterator<IT,T> tmp=*this; iter--; return tmp; }
    inline const rlazyiterator<IT,T>& operator--() // prefix
      { --iter; return *this; }

    inline const rlazyiterator<IT,T>& operator+=(int n)
      { iter+=n; return *this; }
    inline const rlazyiterator<IT,T>& operator-=(int n)
      { iter-=n; return *this; }

//      template <class ITF, class TF> friend const rlazyiterator<ITF,TF> 
//      operator+(const rlazyiterator<ITF,TF>& it, int n)
//        { return rlazyiterator<ITF,TF>(it.iter + n,it.lazyptr); }
//      template <class ITF, class TF> friend const rlazyiterator<ITF,TF> 
//      operator+(int n, const rlazyiterator<ITF,TF>& it)
//        { return rlazyiterator<ITF,TF>(n + it.iter,it.lazyptr); }
//      template <class ITF, class TF> friend const rlazyiterator<ITF,TF> 
//      operator-(const rlazyiterator<ITF,TF>& it, int n)
//        { return rlazyiterator<ITF,TF>(it.iter - n,it.lazyptr); }
//      template <class ITF, class TF> friend const rlazyiterator<ITF,TF> 
//      operator-(int n, const rlazyiterator<ITF,TF>& it)
//        { return rlazyiterator<ITF,TF>(n - it.iter,it.lazyptr); }


    inline bool operator==(const rlazyiterator<IT,T>& it) const
       { return iter == it.iter; }
    inline bool operator!=(const rlazyiterator<IT,T>& it) const 
       { return iter != it.iter; }
    inline bool operator<(const rlazyiterator<IT,T>& it) const
       { return iter < it.iter; }
    inline bool operator>(const rlazyiterator<IT,T>& it) const
       { return iter > it.iter; }
    inline bool operator<=(const rlazyiterator<IT,T>& it) const
       { return iter <= it.iter; }
    inline bool operator>=(const rlazyiterator<IT,T>& it) const
       { return iter >= it.iter; }

    inline const rlazyiterator<IT,T>& 
      operator=(const rlazyiterator<IT,T>& source) 
      { iter=source.iter; lazyptr = source.lazyptr; return *this; }

    inline T& operator*() const 
        { lazyptr->set_whole_cache_validity(false); return *iter;}
    inline T& operator[](int n) const { return *(this + n); }
  };


  //---------------------------------------------------------------------//

  // Constant, Random Access Iterator

  // Use normal constant iterator


  //---------------------------------------------------------------------//

  // Mutable, Bidirectional Iterator

  template <class IT, class T>
  class bilazyiterator {
   private:
    IT iter;
    lazymanager* lazyptr;
  public:
    bilazyiterator() : lazyptr(0) { }
    bilazyiterator(const bilazyiterator<IT,T>& source) : 
      iter(source.iter) , lazyptr(source.lazyptr) { }
    bilazyiterator(const IT& sourceiter, lazymanager* lazyp) :
      iter(sourceiter) , lazyptr(lazyp) { }
    ~bilazyiterator() { }  // do nothing

    inline const bilazyiterator<IT,T> operator++(int) 
      { bilazyiterator<IT,T> tmp=*this; iter++; return tmp; }
    inline const bilazyiterator<IT,T>& operator++() // prefix
      { ++iter; return *this; }
    inline const bilazyiterator<IT,T> operator--(int) 
      { bilazyiterator<IT,T> tmp=*this; iter--; return tmp; }
    inline const bilazyiterator<IT,T>& operator--() // prefix
      { --iter; return *this; }

    inline bool operator==(const bilazyiterator<IT,T>& it) const
       { return iter == it.iter; }
    inline bool operator!=(const bilazyiterator<IT,T>& it) const 
       { return iter != it.iter; }

    inline const bilazyiterator<IT,T>& 
      operator=(const bilazyiterator<IT,T>& source) 
      { iter=source.iter; lazyptr = source.lazyptr; return *this; }

    inline T& operator*() const 
        { lazyptr->set_whole_cache_validity(false); return *iter;}
  };


  //---------------------------------------------------------------------//

  // Constant, Bidirectional Iterator

  // Use normal constant iterator

  //---------------------------------------------------------------------//

  // Mutable, Forward Iterator

  template <class IT, class T>
  class flazyiterator {
  private:
    IT iter;
    lazymanager* lazyptr;
  public:
    flazyiterator() : lazyptr(0) { }
    flazyiterator(const flazyiterator<IT,T>& source) : 
      iter(source.iter) , lazyptr(source.lazyptr) { }
    flazyiterator(const IT& sourceiter, lazymanager* lazyp) :
      iter(sourceiter) , lazyptr(lazyp) { }
    ~flazyiterator() { }  // do nothing

    inline const flazyiterator<IT,T> operator++(int) 
      { flazyiterator<IT,T> tmp=*this; iter++; return tmp; }
    inline const flazyiterator<IT,T>& operator++() // prefix
      { ++iter; return *this; }

    inline bool operator==(const flazyiterator<IT,T>& it) const
       { return iter == it.iter; }
    inline bool operator!=(const flazyiterator<IT,T>& it) const 
       { return iter != it.iter; }

    inline const flazyiterator<IT,T>& 
      operator=(const flazyiterator<IT,T>& source) 
      { iter=source.iter; lazyptr = source.lazyptr; return *this; }

    inline T& operator*() const 
        { lazyptr->set_whole_cache_validity(false); return *iter;}
  };


  //---------------------------------------------------------------------//

  // Constant, Forward Iterator

  // Use normal constant iterator

  //---------------------------------------------------------------------//

}  // end namespace

#endif

