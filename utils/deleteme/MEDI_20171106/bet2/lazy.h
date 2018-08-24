/*  Lazy evaluation support

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

#if !defined(__lazy_h)
#define __lazy_h

#include <iostream>
#include <map>
#include <cstdlib>

#ifndef NO_NAMESPACE
using namespace std;
 namespace LAZY {
#endif


typedef std::map<unsigned int,bool,less<unsigned int> > mapclass;
typedef std::map<unsigned int,bool,less<unsigned int> >::iterator mapiterator;


//-------------------------------------------------------------------------//

template <class T, class S>
class lazy {
private:
  mutable T storedval;
  unsigned int tag;
  const S *iptr;
  T (*calc_fn)(const S &); 

private:
  const T& value() const;
  T calculate_val() const { return (*calc_fn)(*iptr); }

public:
  lazy() { tag = 0; }
  void init(const S *, T (*fnptr)(const S &));
  void copy(const lazy &, const S *);

  const T& force_recalculation() const;

  const T& operator() () const { return this->value(); }
};


//-------------------------------------------------------------------------//

class lazymanager {
  template <class T, class S>  friend class lazy;
private:
  mutable bool validflag;
  mutable mapclass validcache;
  mutable unsigned int tagnum;

private:
  unsigned int getnewtag() const { return tagnum++; }

  const bool is_whole_cache_valid() const 
    { return validflag; }

  const bool is_cache_entry_valid(const unsigned int tag) const 
    { return validcache[tag]; }
  void set_cache_entry_validity(const unsigned int tag, const bool newflag) const 
    { validcache[tag] = newflag; }

  void invalidate_whole_cache() const;

public:
  lazymanager(); 
  void copylazymanager(const lazymanager &);
  void set_whole_cache_validity(const bool newflag) const 
    { validflag = newflag; }
};


//-------------------------------------------------------------------------//

// Body of lazy member functions (put here as cannot simply separate
//   templated definitions into seperate source files if building a library)



template <class T, class S>
const T& lazy<T,S>::value() const 
  {
    if ( (iptr == 0) || (tag==0) ) {
      cerr << "Error: uninitialized lazy evaluation class" << endl;
      exit(-1);
    }
    if (! iptr->is_whole_cache_valid() ) {
      iptr->invalidate_whole_cache();
      iptr->set_whole_cache_validity(true);
    }
    if (! iptr->is_cache_entry_valid(tag)) {
      //cerr << "Calculating value" << endl;
      storedval = calculate_val();
      iptr->set_cache_entry_validity(tag,true);
    }
    return storedval; 
  }

 
template <class T, class S>
const T& lazy<T,S>::force_recalculation() const 
  {
    if ( (iptr == 0) || (tag==0) ) {
      cerr << "Error: uninitialized lazy evaluation class" << endl;
      exit(-1);
    }
    // still process the whole cache vailidity so that this calculation
    //  can get cached correctly
    if (! iptr->is_whole_cache_valid() ) {
      iptr->invalidate_whole_cache();
      iptr->set_whole_cache_validity(true);
    }

    storedval = calculate_val();
    iptr->set_cache_entry_validity(tag,true);

    return storedval; 
  }
 

template <class T, class S>
void lazy<T,S>::init(const S *ip, T (*fnptr)(const S &)) 
  { 
    iptr = ip; 
    calc_fn = fnptr; 
    tag = iptr->getnewtag(); 
    iptr->set_cache_entry_validity(tag,false); 
  }
 

template <class T, class S>
void lazy<T,S>::copy(const lazy &source, const S *ip) {
  storedval = source.storedval;
  tag = source.tag;
  calc_fn = source.calc_fn; 
  // Do NOT copy the same parent class pointer 
  //   (allows parent class to be copied correctly)
  iptr = ip;
}

 
#ifndef NO_NAMESPACE
 }
#endif

#endif


