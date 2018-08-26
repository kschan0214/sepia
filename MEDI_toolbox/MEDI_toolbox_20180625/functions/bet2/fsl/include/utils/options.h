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


#if !defined(OPTIONS_H)
#define OPTIONS_H

#include <stdexcept>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <iterator>
#include <cstdlib>

#define POSIX_SOURCE 1

using namespace std;

namespace Utilities {

  bool string_to_T(bool &b, const std::string& s); 
  bool string_to_T(std::string& d, const std::string& s); 
  bool string_to_T(int& i, const std::string& s); 
  bool string_to_T(float& v, const std::string& s); 
  bool string_to_T(std::vector<int>& vi, const std::string& s);
  bool string_to_T(std::vector<float>& vi, const std::string& s);
  bool string_to_T(std::vector<std::string>& vi, const std::string& s);

  typedef enum argflag { 
    no_argument = 0, requires_argument, optional_argument, requires_2_arguments, 
                     requires_3_arguments, requires_4_arguments, requires_5_arguments
  } ArgFlag;

}

namespace Utilities {


  template<class T,class U>
  std::ostream& operator<<(std::ostream& os, const std::pair<T, U>& o)
  {
    return os << o.first << "," << o.second;
  }
  
  template<class T>
  std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
  {
    std::ostringstream oss;
    std::copy(v.begin(), v.end(), std::ostream_iterator<T>(oss,","));
    
    std::string s(oss.str());
  
    os << s.substr(0, s.length() - 1); // Strip off the final ","
    
    return os;
  }

  /**
     Throw this exception if an error occured inside the Options package.
   */
  class X_OptionError: public std::exception {
  public:
    X_OptionError() throw() {}
    X_OptionError(const std::string& o) throw() :
      m_option(o), m_explanation("unknown error") {}
    X_OptionError(const std::string& o, const std::string& e) throw() :
      m_option(o), m_explanation(e) {}

    virtual const char * what() const throw() {
      std::string str(std::string(m_option + ": " + m_explanation + "!"));
      return str.c_str();
    }

    ~X_OptionError() throw() {}

  private:
    std::string m_option;
    std::string m_explanation;
  };

  /**
     Provides behaviour common to all option types. Actual options are
     declared using the templated Option class. The 
     OptionParser class can be used to parse command lines.
     @see Option
     @see OptionParser
     @author Dave Flitney
     @version 1.0b, Nov., 2000.
  */
  class BaseOption {
  public:
    /**
       @param k comma seperated list of key aliases
       @param ht the help text to be printed for this option
       @param c if true then this option is compulsory
       @param f one of no_argument, requires_argument, optional_argument
       to indicate what arguments should be supplied
    */
    BaseOption(const std::string& k, const std::string& ht, bool c, ArgFlag f): 
      key_(k), help_text_(ht), arg_flag_(f), 
      unset_(true), compulsory_(c), visible_(true) {}
    /**
       @param k comma seperated list of key aliases
       @param ht the help text to be printed for this option
       @param c if true then this option is compulsory
       @param f one of no_argument, requires_argument, optional_argument
       @param v true or false -- display the option in usage
       to indicate what arguments should be supplied
    */
    BaseOption(const std::string& k, const std::string& ht, bool c, ArgFlag f, bool v): 
      key_(k), help_text_(ht), arg_flag_(f), 
      unset_(true), compulsory_(c), visible_(v) {}

    /**
       @return true if the option is compulsory
    */
    bool compulsory() const { return compulsory_; }
    /**
       @return true if the option should be visible
    */
    bool visible() const { return visible_; }
    /**
       @param v true if the option should be visible
    */
    void visible(bool v) { visible_ = v; }
    /**
       @return true if the option requires an argument
    */
    bool required() const { return ( (arg_flag_ == requires_argument) || 
			             (arg_flag_ == requires_2_arguments) ||
			             (arg_flag_ == requires_3_arguments) ||
			             (arg_flag_ == requires_4_arguments) ||
			             (arg_flag_ == requires_5_arguments) ); }
    /**
       @return the number of required arguments 
    */
    int nrequired() const { 
      if ((arg_flag_ == requires_argument) || 
	  (arg_flag_ == optional_argument)) { return 1; } 
      else if (arg_flag_ == requires_2_arguments) { return 2; }
      else if (arg_flag_ == requires_3_arguments) { return 3; }
      else if (arg_flag_ == requires_4_arguments) { return 4; }
      else if (arg_flag_ == requires_5_arguments) { return 5; }
      else return 0;
    }
    /**
       @return true if the option has an optional argument
    */
    bool optional() const { return arg_flag_ == optional_argument; }
    /**
       @return true if the option has an argument at all
    */
    bool has_arg() const { return arg_flag_ != no_argument; }
    /**
       @return true if the option has been set
    */
    bool set() const { return !unset_; }
    /**
       @return true if the option remains unset
    */
    bool unset() const { return unset_; }
    /*
      @param arg A command line argument to be compared against
      the list of possible keys for this option.
      @return True if a match is found.
    */
    bool matches(const std::string& arg);
    /*
      @return This options key string.
    */
    const std::string& key() const { return key_; }

    virtual std::string value_string() const = 0;
    virtual std::string config_key() const = 0;

    /*
      @return This options help text.
    */
    const std::string& help_text() const { return help_text_; }

    /*
      @param Sets the value for this option. Is overridden in the type
      specific template class Option.
    */
    virtual bool set_value(const std::string& vs) = 0;
    virtual bool set_value(const std::string& vs, char *argv[], int valpos, int argc) = 0;

    // For use with optional switch arguments...
    bool use_default_value() {
      unset_ = false;
      return true;
    }

    void usage(std::ostream& os) const;

    virtual std::ostream& print(std::ostream& os) const = 0;

    virtual ~BaseOption() {}

  private:
    std::string key_, help_text_;
    ArgFlag arg_flag_;

  protected:
    /*
      @return This options short-form key (if any)
    */
    const std::string short_form() const;
    /*
      @return This options long-form key (if any)
    */
    const std::string long_form() const;

    bool unset_, compulsory_, visible_;
  };


  /**
     Template class adding type specific behaviour to BaseOption. Define
     one of these per program supported option.
     @author Dave Flitney
     @version 1.0b, Nov., 2000.
     @see BaseOption
  */
  template<class T> class Option: public BaseOption {
  public:
    /** 
	@param k Comma seperated list of key aliases
	@param v Default value for this option
	@param ht Help text to be printed when outputting usage
	@param c If true then this option is compulsory
	@param f This options argument requirements
    */
    Option(const std::string& k, const T& v, const std::string& ht,
	   bool c, ArgFlag f = no_argument): 
      BaseOption(k, ht, c, f), default_(v), value_(v) {}
  /** 
	@param k Comma seperated list of key aliases
	@param v Default value for this option
	@param ht Help text to be printed when outputting usage
	@param c If true then this option is compulsory
	@param f This options argument requirements
    */
    Option(const std::string& k, const T& v, const std::string& ht,
	   bool c, ArgFlag f, bool vis): 
      BaseOption(k, ht, c, f, vis), default_(v), value_(v) {}


    /** 
	@param s The value string which needs to be parsed to set
	this options value. The overloaded function string_to_T must be defined
	for type T.

	@return true if the value actually got set
    */
    bool set_value(const std::string& s)
      { 
	if(string_to_T(value_, s))
	  unset_ = false;
	return !unset_;
      }

    // and a version for multiple options...
    bool set_value(const std::string& vs, char* argv[], int valpos, int argc) {

      if (nrequired()<=0) { /* error */ return false; }
      if (nrequired()==1) {
	// first and only argument
	if(string_to_T(value_, vs)) {
	  unset_ = false;
	  return !unset_;
	}
      } else {
	// Multiple argument case
	T tmpval;
	valuevec_size_ = 0;
	std::string vstmp;
	for (int nv=0; nv<nrequired(); nv++) {
	  if (valpos+nv<argc) vstmp=std::string(argv[valpos+nv]); else vstmp=std::string();
	  unset_ = !string_to_T(tmpval,vstmp);
	  valuevec_[nv] = tmpval;
	  valuevec_size_ = nv + 1;
	}
	return !unset_;
      }
      return !unset_;
    }

    std::string config_key() const
    {
      std::string key(long_form());
      if( key != "" )
	key = key + "=";
      else
	key = short_form() + " ";

      return key;
    }

    std::string value_string() const
    {
      std::ostringstream os;
      os << value();

      return os.str();
    }

    /** 
	@param v The new value.
	Set the options type specific value. You should rarely need this
	method. Use of it means you either haven't planned your options well
	or you have some very complex option rules. You should probably issue
	a warning in your user interface if you call this as you're possibly
	overriding his/her prefered value.
    */
    bool set_T(const T& v) { 
      value_ = v;
      unset_ = false;
      return !unset_;
    }

    /**
       @return The options value.
    */
    const T& value(int n=0) const { 
      if (nrequired()<=1) { return value_; }
      else if ((unsigned int) n<=valuevec_size_) {
	return valuevec_[n]; 
      } else {
	throw X_OptionError("Request for more arguments in value() than exist");
      }
    }
    /** 
	@return The options default value.
    */
    const T& default_value() const { return default_; }

    virtual ostream& print(ostream& os) const {
      os << "# " << help_text() << std::endl 
	 << config_key() << value_string();
      
      return os;
    }

    virtual ~Option() {}

  protected:
    void display() const;

  private:
    Option() {}
  
    T default_, value_;
    T valuevec_[5];
    unsigned int valuevec_size_;
  };

  /*template<> string Option<bool>::config_key() const;
  template<> string Option<bool>::value_string() const; 
  template<> bool Option<bool>::set_value(const string& s);
  template<> std::ostream& Option<bool>::print(std::ostream& s) const;*/
  //  std::ostream& operator<<(std::ostream& os, const Option<bool>& o);
  std::ostream& operator<<(std::ostream& os, const BaseOption& o);

  template<class T> class HiddenOption: public Option<T>
  {
  public:
    HiddenOption(const std::string& k, const T& v, const std::string& ht,
		 bool c, ArgFlag f = no_argument): 
      Option<T>(k, v, ht, c, f, false) {}
  };
  
  template<class T> class FmribOption: public Option<T>
  {
  public:
    FmribOption(const std::string& k, const T& v, const std::string& ht,
		 bool c, ArgFlag f = no_argument): 
      Option<T>(k, v, ht, c, f) { if( getenv("FSLINFMRIB") ){ Option<T>::visible_ = true; } else { Option<T>::visible_ = false; } }
  };
  
  /**
     A class for parsing command line arguments into Option objects. The 
     following example demonstrates the principle:
     @author Dave Flitney
     @version 1.0b, Nov., 2000.
     @see BaseOption
     @see Option
     <pre>
#include "options.h"

// $Id: options.h,v 1.33 2009/04/24 15:34:42 flitney Exp $ 

using namespace Utilities;

Option<bool> verbose(string("-V,--verbose"), false, 
		     string("switch on diagnostic messages"), 
		     false, BaseOption::no_argument);
Option<bool> help(string("-h,--help"), false,
		  string("display this message"),
		  false, BaseOption::no_argument);
Option<float> dof(string("-d,--dof"), 100.0,
		  string("number of degrees of freedom"),
		  true, BaseOption::requires_argument);
Option<string> mask(string("-m,--mask"), string("mask"),
		    string("brain mask volume"),
		    true, BaseOption::requires_argument);
Option<string> resid(string("-r,--res"), string("res4d"),
		     string("4d `residual-of-fit' image"),
		     true, BaseOption::requires_argument);
Option<int> segid(string("-s,--shared-seg-id"), -1,
		  "shared memory segment ID",
		  false, BaseOption::requires_argument);

string title = 
"opttst (Version 1.0)\n\n\
Copyright(c) 2000, University of Oxford\n\
Dave Flitney";

string examples =
"opttst --dof=<number> --mask=<filename> --res=<filename>\n\
opttst -d <number> -m <filename> -r <filename>\n\
opttst --verbose\n";

int main(unsigned int argc, char **argv) {

  OptionParser options(title, examples);

  try {

    options.add(verbose);
    options.add(help);
    options.add(segid);
    options.add(dof);
    options.add(mask);
    options.add(resid);

    for(unsigned int a = options.parse_command_line(argc, argv); 
	a < argc; a++)
      cout << argv[a] << endl;
  
    if(help.value() || 
       !options.check_compulsory_arguments())
      options.usage();

    if(verbose.value()) {
      cout << "verbose = " << verbose.value() << endl;
      cout << "help = " << help.value() << endl;
      cout << "segid = " << segid.value() << endl;
      cout << "dof = " << dof.value() << endl;
      cout << "mask = " << mask.value() << endl;
      cout << "resid = " << resid.value() << endl;
    }

  } catch(X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
  } catch(std::exception &e) {
    cerr << e.what() << endl;
  }    
}
     </pre>
  */
  class OptionParser {
  public:

  OptionParser(const std::string& p, const std::string& e): progname_(p), example_(e), overWriteMode_(Allow) {}

    /**
       @param o An option to be added to the parser
    */
    void add(BaseOption& o) { options_.push_back(&o); }

    void describe_options();
    void brief_usage();
    void usage();
    /**
       @param verbose If set then this method will carp about any option which
       is marked as compulsory but hasn't been set
       @return true if all compulsory arguments have been set and false otherwise
    */
    bool check_compulsory_arguments(bool verbose=false);
    /**
       The parameters, argc and argv, should normally be those passed to
       main via the command line shell.
       @param argc The argument count.
       @param argv The vector of argument strings.
    */
    unsigned int parse_command_line(unsigned int argc, char **argv, int skip=0);

    /**
       @param filename The config file name.
    */
    unsigned int parse_config_file(const std::string& filename);

    virtual ~OptionParser() {}

  protected:
    /**
    */
    OptionParser() {}

  private:
    enum OverwriteMode {Allow=0, ThrowException, Ignore};

    /**
       @param optstr A string which should match one of the option strings
       registered with the add method.
       @return Pointer to the matching option or NULL if a match wasn't found.
    */
    BaseOption* find_matching_option(const std::string& optstr);
    /**
       @param optstr A string which should match one of the option strings
       registered with the add method.
       @param valstr A string which can be used to set the options value
       if applicable.
       @return true on success.
    */
    unsigned int parse_option(const std::string& optstr, const std::string& valstr, 
			      char *argv[], int valpos, int argc)
      throw(X_OptionError);

    /**
       @param str A string of the form --option[=value].
       @return true on success.
    */
    unsigned int parse_long_option(const std::string& str);

    std::string progname_, example_;

    typedef std::vector<BaseOption *> Options;
    Options options_;
    OverwriteMode overWriteMode_;

    friend std::ostream& operator<<(std::ostream& os, const OptionParser p);
  };

}

// std::ostream& operator<<(std::ostream& os, const std::pair<float, float>& o)
// {
//   return os << o.first << "," << o.second;
// }

#endif
