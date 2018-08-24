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
#include <fstream>

namespace Utilities {

  using namespace std;

  BaseOption * OptionParser::find_matching_option(const string& optstr)
  {
    for(Options::iterator o = options_.begin();
	o != options_.end();
	++o)
      if((*o)->matches(optstr))
	return *o;

    return 0;
  }

  unsigned int OptionParser::parse_option(const string& optstr, const string& valstr,
					  char *argv[], int valpos, int argc)
    throw(X_OptionError)
  {
    BaseOption * theOption = 0;

    if((theOption = find_matching_option(optstr)) == 0)
      throw X_OptionError(optstr, "Option doesn't exist");

    if(theOption->unset() || (overWriteMode_==Allow)) 
      {
	if(theOption->has_arg()) {
	  if(valstr.length() > 0) {
	    if(theOption->set_value(valstr,argv,valpos,argc))
	      return 1 + theOption->nrequired();
	    else {
	      string errstr = "Couldn't set_value! valstr=\"" + valstr;
	      for (int nn=valpos+1; nn<=valpos + theOption->nrequired(); nn++) {
		if (nn<argc)  errstr += " " + string(argv[nn]);
	      }
	      throw X_OptionError(optstr, errstr + "\""); 
	    }
	  } else if(!theOption->optional()) {
	    throw X_OptionError(optstr, "Missing non-optional argument");
	  }
	}
	if(theOption->optional()) 
	  theOption->use_default_value();
	else
	  theOption->set_value(string());
	return 1;
      } 
    else 
      {
	if( overWriteMode_!= Ignore)
	  throw X_OptionError(optstr, "Option already set");
	else
	  return 1;
      }

    throw X_OptionError(optstr);
    return 0;
  }


  unsigned int OptionParser::parse_long_option(const string& str)
  {
    string optstr(str);
    string valstr;

    string::size_type pos = 0;
    if((pos = str.find("=", 0)) != string::npos) {
      optstr = str.substr(0, pos);
      valstr = str.substr(pos + 1, str.length() - pos + 1);
    }

    parse_option(optstr, valstr, 0,0,0);

    return 1;
  }

  unsigned int OptionParser::parse_config_file(const string& filename)
  {
    ifstream cf(filename.c_str());

    if(cf.fail())
      throw X_OptionError(filename, "Couldn't open the file");
    
    OverwriteMode oldMode=overWriteMode_;
    overWriteMode_=Ignore;

    string optstr; char buffer[1024];
    while (cf >> optstr) {
      if(optstr[0] == '#')
	cf.getline(buffer, 1024);	     // Read and discard the rest of this line
      else if(optstr.substr(0,2) == "--")
	parse_long_option(optstr); // Parse a long option
      else {
	cf.getline(buffer, 1024);
	parse_option(optstr, string(buffer), 0, 0, 0);
      }
    }
    overWriteMode_=oldMode;
    return 1;
  }
 
  unsigned int OptionParser::parse_command_line(unsigned int argc, 
						char **argv, int skip) 
  {
    unsigned int optpos = 1 + skip;
    unsigned int valpos = 1 + skip;

    while(optpos < argc) {

      unsigned int increments = 0;
      
      string optstr(argv[optpos]), valstr;

      if(optstr[0] != '-')	// End of parsable options
	break;

      if(optstr[1] == '-') {	// Parse a long opt

	increments = parse_long_option(optstr);
	optpos += increments;

      } else {

	valpos = optpos + 1;

	for(unsigned int i = 1; i < optstr.length(); ++i)
	  {
	    string suboptstr = "-" + optstr.substr(i, 1);
	    
	    if (valpos<argc) valstr=string(argv[valpos]); else valstr=string();
	    increments = parse_option(suboptstr, valstr, argv, valpos, argc);
	    
	    valpos += increments - 1;
	  }
	
	optpos = valpos;
      }
    } 
    return optpos;		// User should process any remaining args
  }

  std::ostream& operator<<(std::ostream& os, const OptionParser p) 
  {
    for(OptionParser::Options::const_iterator o = p.options_.begin();
	o != p.options_.end(); ++o)
      os << *(*o) << std::endl;

    return os;
  }
}
