/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */

#include "buildLIZandCommLists.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

#include "Real.hpp"
#include "Neighbors.hpp"

/** sort_atoms function, David M. Rogers
 *
 * Example neighbor-list construction.
 *
 *  const int nx=4, ny=4, nz=4;
 *  const Real Rc = 4.0;
 *  Geom geo(10.0, 12.0, 11.0, nx, ny, nz);
 *  auto cell = sort_atoms(crystal, geo);
 *  NeighborCells nbr(geo, Rc);
 *
 *  for(int ci=0; ci<geo.cells(); ci++) {
 *      int ai, aj, ak;
 *      geo.decodeBin(ci, ai, aj, ak);
 *      for(const LatticePt &pt : nbr) {
 *          int cj = geo.calcBin((pt.i+ai)%geo.n[0],
 *                               (pt.j+aj)%geo.n[1],
 *                               (pt.k+ak)%geo.n[2]);
 *          for(const auto i : cell[ci]) {
 *              for(const auto j : cell[cj]) {
 *                  if(distance < Rc) {
 *                      nbr_list[i].push_back(j);
 *                  }
 *              }
 *          }
 *      }
 *  }
 */
std::vector<std::vector<int>> sort_atoms(CrystalParameters &crystal, const Geom &geo) {
    std::vector<std::vector<int>> cell( geo.cells() );

    for(int a=0; a<crystal.num_atoms; a++) {
        Real x = crystal.position(0, a);
        Real y = crystal.position(1, a);
        Real z = crystal.position(2, a);
        auto ci = geo.calcBinF(x, y, z);
        cell[ci].push_back(a);
    }
    return cell;
}

// Relative brick coordinates shifts x,y,z relative
// to their closest brick origin. - written by David M. Rogers
inline void relative_brick_crd(const Geom &geo, Real &x, Real &y, Real &z) {
    // 1. wrap into range [0,Lz)
    const Real wz = floor(z / geo.L[Geom::ZZ]);
    z -= wz*geo.L[Geom::ZZ];
    y -= wz*geo.L[Geom::ZY];
    x -= wz*geo.L[Geom::ZX];
    const int nz = z/geo.h[Geom::ZZ];

    // 2. wrap into range y_0(nz) + [0,Ly)
    const Real y0 = nz*geo.h[Geom::ZY];
    const auto wy = floor((y-y0) / geo.L[Geom::YY]);
    y -= wy*geo.L[Geom::YY];
    x -= wy*geo.L[Geom::YX];
    const int ny = (y - y0)/geo.h[Geom::YY];

    // 3. wrap into range x_0(ny,nz) + [0,Lx)
    const Real x0 = nz*geo.h[Geom::ZX] + ny*geo.h[Geom::ZY];
    x -= geo.L[Geom::XX]*floor((x - x0)/geo.L[Geom::XX]);
    const int nx = (x - x0)/geo.h[Geom::XX];

    // make x, y, z positions relative to box base crd.
    z -= nz*geo.h[Geom::ZZ];
    y -= ny*geo.h[Geom::YY] + y0;
    x -= nx*geo.h[Geom::XX] + x0;
}

/*
void setupVPboundaries(CrystalParameters &crystal, int idx, AtomData &atom,
		       const Geom &geo, const std::vector<std::vector<int>> &cell)
{
    const Real rtol = 1.0e-8;
    const Real rcirclu = crystal.types[crystal.type[idx]].rLIZ;
    const Real rcirclusqr = rcirclu*rcirclu;
    NeighborCells nbr(geo, rcirclu);
    int nrsclu = 0;


    Real x = crystal.position(0,idx);
    Real y = crystal.position(1,idx);
    Real z = crystal.position(2,idx);
    auto bin = geo.calcBinF(x, y, z); // bin of atom idx
    int ai, aj, ak; // bin's LatticePt
    geo.decodeBin(bin, ai, aj, ak);
    ai += geo.n[0]; aj += geo.n[1]; ak += geo.n[2];

    relative_brick_crd(geo, x, y, z); // wrap relative to bin

    for(const LatticePt &pt : nbr) {
        int cj = geo.calcBin((ai+pt.i)%geo.n[0],
                             (aj+pt.j)%geo.n[1],
                             (ak+pt.k)%geo.n[2]);

        Real offset[3];
        geo.offset(pt, offset); // calculate offset to far cell base pt.
        offset[0] -= x; // make relative to ptcle at idx
        offset[1] -= y;
        offset[2] -= z;

        for(const int j : cell[cj]) {
            Real dx = crystal.position(0,j);
            Real dy = crystal.position(1,j);
            Real dz = crystal.position(2,j);
            relative_brick_crd(geo, dx, dy, dz);
            dx += offset[0]; // add relative bin shifts
            dy += offset[1];
            dz += offset[2];

            const Real atdistsqr = dx*dx + dy*dy + dz*dz;
            if(atdistsqr <= rcirclusqr) {
                LIZ[nrsclu].idx = j; // far atom index
                LIZ[nrsclu].p1 = dx; LIZ[nrsclu].p2 = dy; LIZ[nrsclu].p3 = dz;
                LIZ[nrsclu++].dSqr = atdistsqr;
            }
        }
    }
    return nrsclu;
}
*/

/* O(1) Cell-based neighbor list construction
 * written by David M. Rogers
 */
int buildLIZ(CrystalParameters &crystal, int idx, std::vector<LIZInfoType> &LIZ,
             const Geom &geo, const std::vector<std::vector<int>> &cell) {
    const Real rtol = 1.0e-8;
    const Real rcirclu = crystal.types[crystal.type[idx]].rLIZ;
    const Real rcirclusqr = rcirclu*rcirclu;
    NeighborCells nbr(geo, rcirclu);
    int nrsclu = 0;

    if(std::abs(rcirclu) < rtol) { // special case: a one atom LIZ
        LIZ[0].idx = idx; LIZ[0].dSqr = 0.0;
        LIZ[0].p1 = 0.0; LIZ[0].p2 = 0.0; LIZ[0].p3 = 0.0;
        nrsclu = 1;
        return nrsclu;
    }

    Real x = crystal.position(0,idx);
    Real y = crystal.position(1,idx);
    Real z = crystal.position(2,idx);
    auto bin = geo.calcBinF(x, y, z); // bin of atom idx
    int ai, aj, ak; // bin's LatticePt
    geo.decodeBin(bin, ai, aj, ak);
    ai += geo.n[0]; aj += geo.n[1]; ak += geo.n[2];

    relative_brick_crd(geo, x, y, z); // wrap relative to bin

    for(const LatticePt &pt : nbr) {
        int cj = geo.calcBin((ai+pt.i)%geo.n[0],
                             (aj+pt.j)%geo.n[1],
                             (ak+pt.k)%geo.n[2]);

        Real offset[3];
        geo.offset(pt, offset); // calculate offset to far cell base pt.
        offset[0] -= x; // make relative to ptcle at idx
        offset[1] -= y;
        offset[2] -= z;

        for(const int j : cell[cj]) {
            Real dx = crystal.position(0,j);
            Real dy = crystal.position(1,j);
            Real dz = crystal.position(2,j);
            relative_brick_crd(geo, dx, dy, dz);
            dx += offset[0]; // add relative bin shifts
            dy += offset[1];
            dz += offset[2];

            const Real atdistsqr = dx*dx + dy*dy + dz*dz;
            if(atdistsqr <= rcirclusqr) {
                LIZ[nrsclu].idx = j; // far atom index
                LIZ[nrsclu].p1 = dx; LIZ[nrsclu].p2 = dy; LIZ[nrsclu].p3 = dz;
                LIZ[nrsclu++].dSqr = atdistsqr;
            }
        }
    }
    return nrsclu;
}

int buildLIZ(CrystalParameters &crystal, int idx,std::vector<LIZInfoType> &LIZ)
{
  const Real rtol=1.0e-8;
  const int n_max=5;
  int nrsclu=0;
  Real r1,r2,r3,p1,p2,p3,atdistsqr,shift_1,shift_2,shift_3;
  Real rcirclu=crystal.types[crystal.type[idx]].rLIZ;
  Real rcirclusqr=rcirclu*rcirclu;
  int n0,n1,n2,n3;
  int i[(2*n_max+1)*(2*n_max+1)*(2*n_max+1)];
  int j[(2*n_max+1)*(2*n_max+1)*(2*n_max+1)];
  int k[(2*n_max+1)*(2*n_max+1)*(2*n_max+1)];

  r1=std::sqrt(crystal.bravais(0,0)*crystal.bravais(0,0) +
               crystal.bravais(1,0)*crystal.bravais(1,0) +
               crystal.bravais(2,0)*crystal.bravais(2,0));
  r2=std::sqrt(crystal.bravais(0,1)*crystal.bravais(0,1) +
               crystal.bravais(1,1)*crystal.bravais(1,1) +
               crystal.bravais(2,1)*crystal.bravais(2,1));
  r3=std::sqrt(crystal.bravais(0,2)*crystal.bravais(0,2) +
               crystal.bravais(1,2)*crystal.bravais(1,2) +
               crystal.bravais(2,2)*crystal.bravais(2,2));

/*
  n1=std::max(n_max,int(rcirclu/r1+0.9));
  n2=std::max(n_max,int(rcirclu/r2+0.9));
  n3=std::max(n_max,int(rcirclu/r3+0.9));
*/
  n1=std::max(1,int(rcirclu/r1+0.9));
  n2=std::max(1,int(rcirclu/r2+0.9));
  n3=std::max(1,int(rcirclu/r3+0.9));

  if((2*n1+1)*(2*n2+1)*(2*n3+1) > (2*n_max+1)*(2*n_max+1)*(2*n_max+1))
  {
    fprintf(stderr,"FATAL ERROR: buildLIZ: n0 exeeds upper limit for site %d!",idx);
    exit(1);
  }
  n0=0;
  for(int i0=-n1; i0<=n1; i0++)
    for(int j0=-n2; j0<=n2; j0++)
      for(int k0=-n3; k0<=n3; k0++)
      {
        i[n0]=i0; j[n0]=j0; k[n0]=k0; n0++;
      }
  if(std::abs(rcirclu)<rtol) // special case: a one atom LIZ
  {
    LIZ[0].idx=idx; LIZ[0].dSqr=0.0;
    LIZ[0].p1=0.0; LIZ[0].p2=0.0; LIZ[0].p3=0.0;
    nrsclu=1;
    return nrsclu;
  } else {
    nrsclu=0;
    for(int m=0; m<n0; m++)
    {
      shift_1=Real(i[m])*crystal.bravais(0,0) +
        Real(j[m])*crystal.bravais(0,1) +
        Real(k[m])*crystal.bravais(0,2) -
        crystal.position(0,idx);
      shift_2=Real(i[m])*crystal.bravais(1,0) +
        Real(j[m])*crystal.bravais(1,1) +
        Real(k[m])*crystal.bravais(1,2) -
        crystal.position(1,idx);
      shift_3=Real(i[m])*crystal.bravais(2,0) +
        Real(j[m])*crystal.bravais(2,1) +
        Real(k[m])*crystal.bravais(2,2) -
        crystal.position(2,idx);
      for(int n=0; n<crystal.num_atoms; n++)
      {
        p1=shift_1+crystal.position(0,n);
        p2=shift_2+crystal.position(1,n);
        p3=shift_3+crystal.position(2,n);
        atdistsqr=p1*p1+p2*p2+p3*p3;
        if(atdistsqr<=rcirclusqr)
        {
          LIZ[nrsclu].idx=n;
          LIZ[nrsclu].p1=p1; LIZ[nrsclu].p2=p2; LIZ[nrsclu].p3=p3;
          LIZ[nrsclu++].dSqr=atdistsqr;
        }
      }
    }
  }
  return nrsclu;
}

bool nodeIsInList(int val, std::vector<LIZInfoType> &list, int len, CrystalParameters &crystal, int &ret)
{
  for(int i=ret; i<len; i++) if(crystal.types[crystal.type[list[i].idx]].node==val) {ret=i; return true;}
  return false;
}


 
bool nodeLess_NodeIndexInfo(const NodeIdxInfo &x, const NodeIdxInfo &y) {return x.node<y.node;}
bool globalLess_NodeIndexInfo(const NodeIdxInfo &x, const NodeIdxInfo &y) {return x.globalIdx<y.globalIdx;}
bool localLess_NodeIndexInfo(const NodeIdxInfo &x, const NodeIdxInfo &y) {return x.localIdx<y.localIdx;}
bool globalEq_NodeIndexInfo(const NodeIdxInfo &x, const NodeIdxInfo &y) {return x.globalIdx==y.globalIdx;}
bool globalAndNodeEq_NodeIndexInfo(const NodeIdxInfo &x, const NodeIdxInfo &y)
{return (x.globalIdx==y.globalIdx) && (x.node==y.node);}
bool localAndNodeEq_NodeIndexInfo(const NodeIdxInfo &x, const NodeIdxInfo &y)
{return (x.localIdx==y.localIdx) && (x.node==y.node);}
bool dSqrLess_LIZInfoType(const LIZInfoType &x, const LIZInfoType &y) {return x.dSqr<y.dSqr;}

void buildLIZandCommLists(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                          CrystalParameters &crystal, LocalTypeInfo &local)
{
  int numNodes = comm.size;
  std::vector<NodeIdxInfo> toList, fromList;
  std::vector<LIZInfoType> tempLIZ;
  tempLIZ.resize(4096);
  int tempNumLIZ;
  int ret;
  int fromListN,toListN;
  int toCounts[4096];
  int fromCounts[4096];
  int num_store=local.num_local;

  double timeBuildLIZandCommList = MPI_Wtime();
  // Begin cell sorting
  const Real Rc_target = 4.0; // set this to determine nx, ny, nz
  // Real Rc_target = 4.0;
  // for(int i=0; i<local.num_local; i++)
  //   if(Rc_target < local.atom[i].rLIZ)
  //     Rc_target = local.atom[i].rLIZ;
  
  int nx = crystal.bravais(0,0)/Rc_target;
  int ny = crystal.bravais(1,1)/Rc_target;
  int nz = crystal.bravais(2,2)/Rc_target;
  if(nx < 1) nx = 1;
  if(ny < 1) ny = 1;
  if(nz < 1) nz = 1;

  Geom geo(crystal.bravais(0,0), crystal.bravais(1,1), crystal.bravais(2,2),
           nx, ny, nz,
           crystal.bravais(0,1), crystal.bravais(0,2), crystal.bravais(1,2));
  // Note: crystal.bravais(1,0), crystal.bravais(2,0), crystal.bravais(2,1)
  //       must all be zero, or a whole-cell rotation is needed to axis-align a and b
  //       into the x-axis and x-y plane, respectively!
  auto cell = sort_atoms(crystal, geo);
  // End cell sorting
  timeBuildLIZandCommList = MPI_Wtime() - timeBuildLIZandCommList;
  if (lsms.global.iprint >= 0)
  {
    if(lsms.lsmsAlgorithms[LSMSAlgorithms::lizConstruction] == LSMSAlgorithms::lizConstruction_bricks)
      printf("lizConstruction algorithm : bricks\n");
    else printf("lizConstruction algorithm : original\n");

    printf("  time for cell sorting: %lf sec\n",
           timeBuildLIZandCommList);
    fflush(stdout);
  }
  timeBuildLIZandCommList = MPI_Wtime();

  // the first num_local entries in tmatStore contain the local tmats
  local.tmatStoreGlobalIdx.resize(4096);

  for(int i=0; i<local.num_local; i++)
  {
    crystal.types[local.global_id[i]].store_id=i;
    local.tmatStoreGlobalIdx[i]=local.global_id[i];
  }

  for(int i=0; i<4096; i++) {toCounts[i]=fromCounts[i]=0;}

  fromListN=toListN=0;
  fromList.resize(4096*local.num_local);
  toList.resize(4096*local.num_local);
// loop over all sites:
  for(int i=0; i<crystal.num_atoms; i++)
  {
    int type_id=crystal.type[i];
    int node=crystal.types[type_id].node;
    
    
    if(node==comm.rank) // local atom type
    {
      if(lsms.lsmsAlgorithms[LSMSAlgorithms::lizConstruction] == LSMSAlgorithms::lizConstruction_bricks)
      {
        //tempNumLIZ=buildLIZ(crystal,i,tempLIZ);
        tempNumLIZ=buildLIZ(crystal,i,tempLIZ, geo, cell);
      } else {
        tempNumLIZ=buildLIZ(crystal,i,tempLIZ);
      }

// set LIZ
      int local_id=crystal.types[type_id].local_id;
      std::sort(tempLIZ.begin(),tempLIZ.begin()+tempNumLIZ,dSqrLess_LIZInfoType);
      local.atom[local_id].numLIZ=tempNumLIZ;
      local.atom[local_id].LIZGlobalIdx.resize(tempNumLIZ);
      local.atom[local_id].LIZStoreIdx.resize(tempNumLIZ);
      local.atom[local_id].LIZDist.resize(tempNumLIZ);
      local.atom[local_id].LIZlmax.resize(tempNumLIZ);
      local.atom[local_id].LIZPos.resize(3,tempNumLIZ);
      local.atom[local_id].nrmat=0;
      for(int j=0; j<tempNumLIZ; j++)
      {
        local.atom[local_id].LIZGlobalIdx[j]=tempLIZ[j].idx;
        local.atom[local_id].LIZDist[j]=std::sqrt(tempLIZ[j].dSqr);
        local.atom[local_id].LIZPos(0,j)=tempLIZ[j].p1;
        local.atom[local_id].LIZPos(1,j)=tempLIZ[j].p2;
        local.atom[local_id].LIZPos(2,j)=tempLIZ[j].p3;
// calculate the lmax for the various shells
        local.atom[local_id].lmax = crystal.types[type_id].lmax;
        int lkeep=crystal.types[type_id].lmax;
        local.atom[local_id].lmax = lkeep;
        for(int n1=0; n1<4; n1++)
          if(local.atom[local_id].LIZDist[j]>crystal.types[type_id].rsteps[n1]) lkeep--;
        local.atom[local_id].LIZlmax[j]=lkeep;
        local.atom[local_id].nrmat+=(lkeep+1)*(lkeep+1);
// add to commTmatFrom
        if(crystal.types[crystal.type[tempLIZ[j].idx]].node!=comm.rank) // need this from remote node
        {
          fromList[fromListN].node=crystal.types[crystal.type[tempLIZ[j].idx]].node;
          fromList[fromListN].localIdx=crystal.types[crystal.type[tempLIZ[j].idx]].local_id;
          fromList[fromListN].globalIdx=crystal.type[tempLIZ[j].idx];
          fromListN++;
        }
      }
    } else { // non local atom
// before building the LIZ we should first check if it is actually needed
// i.e find min distance between atom and local atoms and compare with max liz radius
      //tempNumLIZ=buildLIZ(crystal,i,tempLIZ);
      if(lsms.lsmsAlgorithms[LSMSAlgorithms::lizConstruction] == LSMSAlgorithms::lizConstruction_bricks)
      {
        tempNumLIZ=buildLIZ(crystal,i,tempLIZ, geo, cell);
      } else {
        tempNumLIZ=buildLIZ(crystal,i,tempLIZ);
      }

      ret=0;
      while(nodeIsInList(comm.rank,tempLIZ,tempNumLIZ,crystal,ret)) // the remote node needs the tmat from our node
      {
        toList[toListN].node=crystal.types[crystal.type[i]].node;          // the node that needs a tmat from us
        toList[toListN].localIdx=crystal.types[crystal.type[tempLIZ[ret].idx]].local_id;  // our local index == tmatStore entry
        toList[toListN].globalIdx=crystal.type[i];                         // the type that needs this tmat
        toListN++; ret++;
      }
    }
  }
  timeBuildLIZandCommList = MPI_Wtime() - timeBuildLIZandCommList;
  if (lsms.global.iprint >= 0)
  {
    printf("  time to build LIZs: %lf sec\n",
           timeBuildLIZandCommList);
    fflush(stdout);
  }

  // build the vpClusterGlobalIdx for the Voronoi polyhedra construction
  for(int i=0; i<local.num_local; i++)
  {
    if(local.atom[i].numLIZ > 50)
    {
      local.atom[i].vpClusterGlobalIdx = local.atom[i].LIZGlobalIdx;
      local.atom[i].vpClusterPos.resize(3,local.atom[i].vpClusterGlobalIdx.size());
      for(int j=0; j<local.atom[i].vpClusterGlobalIdx.size(); j++)
      {
	local.atom[i].vpClusterPos(0, j) = local.atom[i].LIZPos(0, j);
	local.atom[i].vpClusterPos(1, j) = local.atom[i].LIZPos(1, j);
	local.atom[i].vpClusterPos(2, j) = local.atom[i].LIZPos(2, j);
      }
    } else {
      // need to construct vpCluster - for the time set it to size zero to use the old algorithm
      local.atom[i].vpClusterGlobalIdx.clear();
    }
  }    

  
  timeBuildLIZandCommList = MPI_Wtime();

// sort toList and fromList
  std::sort(fromList.begin(),fromList.begin()+fromListN,globalLess_NodeIndexInfo);
  std::sort(toList.begin(),toList.begin()+toListN,localLess_NodeIndexInfo);
  std::stable_sort(toList.begin(),toList.begin()+toListN,nodeLess_NodeIndexInfo);
// remove duplicates (only works on sorted lists!):
// for FROM remove all duplicate atom types
  std::vector<NodeIdxInfo>::iterator it=std::unique(fromList.begin(),fromList.begin()+fromListN,globalEq_NodeIndexInfo);
  fromList.resize(it-fromList.begin());
// for TO remove all identical atoms to the same node
  it=std::unique(toList.begin(),toList.begin()+toListN,localAndNodeEq_NodeIndexInfo);
  toList.resize(it-toList.begin());

//  fromList.resize(fromListN);
//  toList.resize(toListN);
  std::sort(fromList.begin(),fromList.end(),nodeLess_NodeIndexInfo);
  std::sort(toList.begin(),toList.end(),nodeLess_NodeIndexInfo);

// count the nodes in toList and fromList
  if(lsms.global.iprint>0) printf("toList.size()=%zu\n",toList.size());
  int numToNodes=0;
  if(toList.size()>0)
  {
    numToNodes=1;
    int h=toList[0].node;
    for(int i=0; i<toList.size(); i++)
    {
      if(lsms.global.iprint>0) printf("toList[%d].node=%d .localIdx=%d .globalIdx=%d\n",i,toList[i].node,
                              toList[i].localIdx,toList[i].globalIdx);
      // if(h!=toList[i].node) {h=toList[i].node; numToNodes++;} else toCounts[numToNodes-1]++;
      if(h!=toList[i].node) {h=toList[i].node; numToNodes++; toCounts[numToNodes-1]=1;} else toCounts[numToNodes-1]++;
    }
  }

  if(lsms.global.iprint>0) printf("fromList.size()=%zu\n",fromList.size());
  int numFromNodes=0;
  if(fromList.size()>0)
  {
    numFromNodes=1;
    int h=fromList[0].node;
    for(int i=0; i<fromList.size(); i++)
    {
      if(lsms.global.iprint>0) printf("fromList[%d].node=%d .localIdx=%d .globalIdx=%d\n",i,fromList[i].node,
                              fromList[i].localIdx,fromList[i].globalIdx);
      if(h!=fromList[i].node) {h=fromList[i].node; numFromNodes++; fromCounts[numFromNodes-1]=1;} else fromCounts[numFromNodes-1]++;
    }
  }

  comm.numTmatTo=numToNodes;
  comm.tmatTo.resize(numToNodes);
  comm.numTmatFrom=numFromNodes;
  comm.tmatFrom.resize(numFromNodes);

  int k=0;
  for(int i=0; i<numToNodes; i++)
  {
    comm.tmatTo[i].tmatStoreIdx.resize(toCounts[i]);
    comm.tmatTo[i].globalIdx.resize(toCounts[i]);
    comm.tmatTo[i].communicationRequest.resize(toCounts[i]);
    comm.tmatTo[i].remoteNode=toList[k].node;
    comm.tmatTo[i].numTmats=toCounts[i];
    for(int j=0; j<toCounts[i]; j++)
    {
      comm.tmatTo[i].globalIdx[j]=local.global_id[toList[k].localIdx];
      comm.tmatTo[i].tmatStoreIdx[j]=toList[k++].localIdx;
    }
  }

  k=0;
  if(lsms.global.iprint>0) printf("numFromNodes=%d\n",numFromNodes);
  for(int i=0; i<numFromNodes; i++)
  {
    if(lsms.global.iprint>0) printf("fromCounts[%d]=%d\n",i,fromCounts[i]);
    comm.tmatFrom[i].tmatStoreIdx.resize(fromCounts[i]);
    comm.tmatFrom[i].globalIdx.resize(fromCounts[i]);
    comm.tmatFrom[i].communicationRequest.resize(fromCounts[i]);
    comm.tmatFrom[i].remoteNode=fromList[k].node;
    comm.tmatFrom[i].numTmats=fromCounts[i];
    for(int j=0; j<fromCounts[i]; j++)
    {
      int g=fromList[k++].globalIdx;
      comm.tmatFrom[i].globalIdx[j]=g;
      if(crystal.types[g].store_id<0)
      {
        local.tmatStoreGlobalIdx[num_store]=g;
        crystal.types[g].store_id=num_store++;
      }
//      if(comm.rank==0)
//        printf("  i=%d j=%d : g=%d num_store=%d crystal.types[g].store_id=%d\n",i,j,g,num_store,crystal.types[g].store_id);
      comm.tmatFrom[i].tmatStoreIdx[j]=crystal.types[g].store_id;
    }
  }
  local.tmatStoreGlobalIdx.resize(num_store);
  int kkrsz2=2*(crystal.maxlmax+1)*(crystal.maxlmax+1);
  local.blkSizeTmatStore=kkrsz2*kkrsz2;
  local.lDimTmatStore=local.blkSizeTmatStore*lsms.energyContour.groupSize();
  local.tmatStore.resize(local.lDimTmatStore,num_store);

// set the StorIdx for the local atom LIZs
  for(int i=0; i<local.num_local; i++)
    for(int j=0; j<local.atom[i].numLIZ; j++)
      local.atom[i].LIZStoreIdx[j]=crystal.types[crystal.type[local.atom[i].LIZGlobalIdx[j]]].store_id;

  timeBuildLIZandCommList = MPI_Wtime() - timeBuildLIZandCommList;
  if (lsms.global.iprint >= 0)
  { 
    printf("  time to build communication lists: %lf sec\n",
           timeBuildLIZandCommList);
    fflush(stdout);
  }

}
