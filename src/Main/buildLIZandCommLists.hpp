
#ifndef LSMS_BUILDLIZANDCOMMLISTS_HPP
#define LSMS_BUILDLIZANDCOMMLISTS_HPP

#include "Real.hpp"

#include "SystemParameters.hpp"
#include "SingleSite/AtomData.hpp"
#include "Communication/LSMSCommunication.hpp"

class LIZInfoType {
public:
  int idx;
  Real p1, p2, p3;
  Real dSqr;
};

class NodeIdxInfo {
public:
  int node, localIdx, globalIdx;
};

void buildLIZandCommLists(LSMSCommunication &comm,
                          LSMSSystemParameters &lsms,
                          CrystalParameters &crystal,
                          LocalTypeInfo &local);


#endif
