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

#include "options.h"

namespace Utilities {
  
  using namespace std;

  template<> string Option<bool>::config_key() const
  {
    if(set()) {
      string key(long_form());
      if( key == "" )
	key = short_form();
      
      return key;
    } else
      return "";
  } 
  
  template<> string Option<bool>::value_string() const { return ""; }

  template<> bool Option<bool>::set_value(const string& s)
  {
    if(s.length() == 0)
      {
	value_ = !default_;
	unset_=false;
      }
    else if (s == "true")
      {
	value_ = true;
	unset_=false;
      }
    else if (s == "false")
      {
	value_ = false;
	unset_=false;
      }
    return !unset_;
  }

  template<> ostream& Option<bool>::print(ostream& os) const 
  {
    os << "# " << help_text() << endl;
    if(set())
      os << config_key().substr(0, config_key().find("=")); 

    return os;
  }
  
  ostream& operator<<(ostream& os, const BaseOption& o)
  {
    return o.print(os);
  }

  bool string_to_T(bool& b, const string& s) {
    b = false;
    return false;
  }

  bool string_to_T(string& d, const string& s) {
    d = s;
    return true;
  }

  bool string_to_T(int& i, const string& s) {
    char *endptr = 0; const char *str = s.c_str();
    i = strtol(str, &endptr, 0);
    if(*endptr == str[s.length()])
      return true;
    else
      return false;
  }

  bool string_to_T(float& v, const string& s) {
    char *endptr = 0; const char *str = s.c_str();
    v = strtod(str, &endptr);
    if(*endptr == str[s.length()])
      return true;
    else
      return false;
  }

  bool string_to_T(vector<int>& vi, const string& s) {
    string str(s), delin(",");
    if(str.find(":")!=string::npos)
      delin = ":";
    str=str+delin;
    vi.clear();
    while(str.size()) {
      int v = atoi(str.substr(0,str.find(delin)).c_str());
      vi.push_back(v);
      str = str.substr(str.find(delin)+1,str.length()-str.find(delin)-1);
    }
    return true;
  }

  bool string_to_T(vector<float>& vi, const string& s) {
    string str(s), delin(",");
    if(str.find(":")!=string::npos)
      delin = ":";
    str=str+delin;
    vi.clear();
    while(str.size()) {
      float v = atof(str.substr(0,str.find(delin)).c_str());
      vi.push_back(v);
      str = str.substr(str.find(delin)+1,str.length()-str.find(delin)-1);
    }
    return true;
  }

  bool string_to_T(vector<string>& vi, const string& s) {
    string str(s), delin(",");
    if(str.find(":")!=string::npos)
      delin = ":";
    str=str+delin;
    vi.clear();
    while(str.size()) {
      string v = str.substr(0,str.find(delin));
      vi.push_back(v);
      str = str.substr(str.find(delin)+1,str.length()-str.find(delin)-1);
    }
    return true;
  }

//   ostream& operator<<(ostream &os, const BaseOption& o) {
//     string test=o.help_text();
//     if ((test.length()>=1) && (test[0]=='~')) {
//       test[0]=' ';
//       return os << "\t" << o.key() << test;
//     } else {
//       return os << "\t" << o.key() << "\t" << o.help_text();
//     }
//   }

  void BaseOption::usage(ostream& os) const {
    string test(help_text());
     if ((test.length()>=1) && (test[0]=='~')) {
       test[0]=' ';
       os << "\t" << key() << test;
     } else {
       os << "\t" << key() << "\t" << help_text();
     }
   }

  bool is_short_form(const string& s)
  {
    return (s.substr(0,2) != "--");
  }


  /*
    @return first short-form key (if any)
  */
  const string BaseOption::short_form() const
  {
    string::size_type pos(0), np;
    
    while( (np = key_.find(",", pos)) != string::npos ) {
      string candidate(key_.substr(pos, np - pos));
      if( is_short_form(candidate) )
	return candidate;
      else
	pos = np + 1;
    }
    string candidate(key_.substr(pos, np - pos));
    if( is_short_form(candidate) )
      return candidate;
    else
      return "";
  }

  /*
    @return first long-form key (if any)
  */
  const string BaseOption::long_form() const
  {
    string::size_type pos(0), np;
    
    while( (np = key_.find(",", pos)) != string::npos ) {
      string candidate(key_.substr(pos, np - pos));

      if( !is_short_form(candidate) )
	return candidate;
      else
	pos = np + 1;
    }
    string candidate(key_.substr(pos, np - pos));
    if( !is_short_form(candidate) )
      return candidate;
    else
      return "";
  }
}
