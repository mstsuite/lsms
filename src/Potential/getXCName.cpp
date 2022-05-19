/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#include "getXCName.hpp"

#include "Main/SystemParameters.hpp"

#include <string>

bool getXCName(LSMSSystemParameters &lsms, std::string &name)
{
  if(lsms.xcFunctional[0]==0) // built in functionals
  {
    switch(lsms.xcFunctional[1])
    {
      case 1: name="von Barth-Hedin (LSMS_1)"; return true;
      case 2: name="Vosko-Wilk-Nusair (LSMS_1)"; return true;
    }
    name="Illegal Exchange-Correlation Functional (built in)!"; return false;
  } else if(lsms.xcFunctional[0]==1) { // libxc functionals
#ifdef USE_LIBXC
    name = lsms.exch_corr->get_name();
    return true;
#else
    name="Illegal Exchange-Correlation Functional (LSMS not linked to libXC)!"; return false;
#endif
  } else if(lsms.xcFunctional[0]==2) { // new functionals
    switch(lsms.xcFunctional[1])
    { 
      // case 1: name="von Barth-Hedin (LSMS_1)"; return true;
      case 2: name="Vosko-Wilk-Nusair (LSMS_3)"; return true;
    }
    name="Illegal Exchange-Correlation Functional (LSMS_3)!"; return false;
  } else { // unknown functional!!
    name="Illegal Exchange-Correlation Functional!"; return false;
  }
  return false;
}
