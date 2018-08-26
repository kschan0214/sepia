/*  Copyright (C) 2000 University of Oxford  */

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

#if !defined(__complexvolume_h)
#define __complexvolume_h

#include "newimage.h" 


namespace NEWIMAGE {

  class complexpoint;

  class complexref {
  private:
    float *m_real;
    float *m_imag;
  public:
    complexref(float* r, float* i) : m_real(r), m_imag(i){}
    ~complexref(){}
    inline float* re_pointer() const { return m_real; }
    inline float* im_pointer() const { return m_imag; }
    const complexpoint& operator=(const complexpoint& val);
  };

  class complexpoint {

  private:
    float m_real;
    float m_imag;

  public:

    complexpoint(){}
    complexpoint(float r, float i){m_real=r; m_imag=i;}
    complexpoint(float r){m_real=r; m_imag=0;}
    complexpoint(const complexref& source){
      m_real = *(source.re_pointer());
      m_imag = *(source.im_pointer());
    }
    ~complexpoint(){}
    float operator=(const float val);
    complexpoint& operator=(const complexpoint& source);
    complexpoint& operator=(const complexref& source);
    inline float re() const { return m_real; }
    inline float im() const { return m_imag; }

    const complexpoint& operator+=(const complexpoint& val); 
    const complexpoint& operator-=(const complexpoint& val); 
    const complexpoint& operator*=(const complexpoint& val); 
    const complexpoint& operator/=(const complexpoint& val);
    complexpoint operator+(const complexpoint& val) const;
    complexpoint operator-(const complexpoint& val) const;
    complexpoint operator*(const complexpoint& val) const;
    complexpoint operator/(const complexpoint& val) const;

    //ostream& operator<<(ostream& s, const complexpoint& val);

    float abs() const;
    float phase() const;
  };

  class complexvolume {
    
  private:
    volume<float> real;
    volume<float> imag;

  public:
    complexvolume(){}
    complexvolume(int xsize, int ysize, int zsize);
    complexvolume(const complexvolume& source);
    complexvolume(const volume<float>& r, const volume<float>& i);
    complexvolume(const volume<float>& r);
    ~complexvolume();
    float operator=(const float val);
    complexvolume& operator=(const complexvolume& source); 
    void destroy();
    int copyproperties(const complexvolume& source);
    int copydata(const complexvolume& source);
 
    const float& re(int x,int y, int z) const { return real(x,y,z); }
    const float& im(int x,int y, int z) const { return imag(x,y,z); }
    float& re(int x,int y, int z) { return real(x,y,z); }
    float& im(int x,int y, int z) { return imag(x,y,z); }

    inline int xsize() const { return real.xsize(); }
    inline int ysize() const { return real.ysize(); }
    inline int zsize() const { return real.zsize(); }
    inline float xdim() const { return real.xdim(); }
    inline float ydim() const { return real.ydim(); }
    inline float zdim() const { return real.zdim(); }
    void setxdim(float x) { real.setxdim(x); imag.setxdim(x); }
    void setydim(float y) { real.setydim(y); imag.setydim(y); }
    void setzdim(float z) { real.setzdim(z); imag.setzdim(z); }
    void setdims(float x, float y, float z){ setxdim(x); setydim(y); setzdim(z); }
    int nvoxels() const { return real.nvoxels(); }
  
    const complexvolume& operator+=(const complexpoint& val);
    const complexvolume& operator-=(const complexpoint& val);
    const complexvolume& operator*=(const complexpoint& val);
    const complexvolume& operator/=(const complexpoint& val);
    const complexvolume& operator+=(const complexvolume& source); 
    const complexvolume& operator-=(const complexvolume& source); 
    const complexvolume& operator*=(const complexvolume& source); 
    const complexvolume& operator/=(const complexvolume& source); 
    complexvolume operator+(const complexpoint& val) const;
    complexvolume operator-(const complexpoint& val) const;
    complexvolume operator*(const complexpoint& val) const;
    complexvolume operator/(const complexpoint& val) const;
    complexvolume operator+(const complexvolume& vol) const;
    complexvolume operator-(const complexvolume& vol) const;
    complexvolume operator*(const complexvolume& vol) const;
    complexvolume operator/(const complexvolume& vol) const;

    inline const std::vector<int>& limits() const { return real.limits(); }
    inline int limits(int n) const { return real.limits(n); }
    inline int minx() const { return real.minx(); }
    inline int maxx() const { return real.maxx(); }
    inline int miny() const { return real.miny(); }
    inline int maxy() const { return real.maxy(); }
    inline int minz() const { return real.minz(); }
    inline int maxz() const { return real.maxz(); }
    inline const std::vector<int>& ROIlimits() const { return real.ROIlimits(); }
    inline int ROIlimits(int n) const { return real.ROIlimits(n); }
    inline bool usingROI() const { return real.usingROI(); }
    inline void setROIlimits(int x0, int y0, int z0, int x1, int y1, int z1) const
      { real.setROIlimits(x0,y0,z0,x1,y1,z1);
      imag.setROIlimits(x0,y0,z0,x1,y1,z1); }
    inline void setROIlimits(const std::vector<int>& lims) const
      { real.setROIlimits(lims);imag.setROIlimits(lims); }
    inline void activateROI() const
      { real.activateROI(); imag.activateROI(); } 
    inline void deactivateROI() const 
      { real.deactivateROI(); imag.deactivateROI(); }

    volume<float> abs() const;
    volume<float> phase() const;
    volume<float>& re();
    volume<float>& im();
    const volume<float>& re() const;
    const volume<float>& im() const;

    complexref operator()(int x,int y, int z)
      { return(complexref(&real(x,y,z),&imag(x,y,z))); }

    
    complexvolume extract_slice(int slice) const;
    void overwrite_slice(const complexvolume& data,int slice);

  };
}
#endif
