#include <vector>
#include <algorithm>
#include <cstdio>

#include "Real.hpp"
#include "Complex.hpp"
#include "Matrix.hpp"
#include "SingleSite/AtomData.hpp"
#include "initializeAtom.hpp"
#include "Core/CoreStates.hpp"
#include "Misc/integrateOneDim.cpp"
#include "Potential/calculateChargesPotential.hpp"

/* initializeAtom(AtomData &a)
   initialize an atom using the following information in AtomData

    Mesh information: xstart, rmt, jmt, jws
    Atom information: ztotss, zcorss, zsemss, zvalss
    (Note that ztotss=zcorss+zsemss+zvalss)
    lmax to limit core search
*/
class InitialAtomLevels {
public:
  Real energy;
  int n,l;
  std::vector<Real> rho; // radial charge density contribution from an electron in this level
  char type; // 'C' for core lectron, 'S' for semi-core and 'V' for valence electrons
// bool operator<()(InitialAtomLevels const &a, InitialAtomLevels const &b) { return a.energy<b.energy; }
} ;

struct compareInitialAtomLevels {
  bool operator()(InitialAtomLevels const &a, InitialAtomLevels const &b) { return a.energy<b.energy; }
};

/*
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine deepst(nqn,lqn,kqn,en,rv,r,rf,h,z,c,
     >                  nitmax,tol,nws,nlast,iter,iprpts,ipdeq)

c               nqn:     principal quantum number; 
c               lqn:     orbital quantum number;
c               kqn:     kappa quantum number; 
c               en:      energy; 
c               rv:      potential in rydbergs times r;
c               r:       radial log. grid; 
c               rg:      big component; 
c               rf:      small component;
c               h:       exp. step; 
c               z:       atomic number; 
c               nitmax:  number of iterations;
c               tol:     energy tolerance; 
c               nws:     bounding sphere radius index;
c               nlast:   last tabulation point; 
c               c:       speed of light in rydbergs;
c               drg,drf: wavefunctions derivatives times r;
c               gam:     first power for small r expansion; 
c               slp:     slope at the origin;
c               dm:      h/720;
*/
extern "C"
{
  void deepst_check_(int *nqn, int *lqn, int *kqn, Real *en, Real *rv, Real*r,
               Real *rf, Real*h, Real *z, Real *c,
               int *nitmax, Real *tol, int *nws, int *nlast, int *iter,
		     int *iprpts, int *ipdeq, int *found);
}

// approximation for the error function
Real approxErfc(Real x)
{
  const Real p  = 0.47047;
  const Real a1 = 0.3480242;
  const Real a2 = -0.0958798;
  const Real a3 = 0.7478556;

  Real t = 1.0 / (1.0 + p*x);
  return 1.0 - ((a1 + (a2 + a3*t) * t) * t) * std::exp(-(x*x));
}

int calculateAtomLevels(AtomData &a, std::vector<InitialAtomLevels> &atomLevels,
			Matrix<Real> &rhoValence, bool writeOrbitals=false)
{
  // find the lowest zcorss+zsemss levels
  // maximum number of states to test: (lmax+1)*(lmax+2)/2
  // assuming kappa degeneracy, iterating over principal quantum number and l
  // we treat core and semi-core states the same for the initialization,
  // using deepst for both

  int lmaxCore = std::min(a.lmax, 3);   // only consider s,p,d,f electrons (we should never need g or higher)

  int  numAtomLevels = ((lmaxCore+1) * (lmaxCore+2)) / 2;
  int  kappa         = 0;
  Real energy        = 0.0;
  Real cLight        = 2.0 * 137.0359895;
  int  nitmax        = 100;                // maximum number of iterations
  Real tol           = 1.0e-8;             // energy tolerance
  int  last          = a.r_mesh.size();
  int  ipdeq         = 5;
  int  iter          = 0;
  int found;

  std::vector<Real> rf(a.r_mesh.size()+1);
  std::vector<Real> rfNormalized(a.r_mesh.size()+1);
  std::vector<Real> unnormalizedDensity(a.r_mesh.size()+1);

  int nlast;
  
  int nl = 0;
  
  for(int n=1; n<=lmaxCore+1; n++)
  {
    for (int l=0; l<n; l++)
    {
      kappa  = -l - 1;        //only use kappa -l-1 this will work for all l (including l=0)
      energy = -(a.ztotss) * (a.ztotss) / Real(n*n);

      if(writeOrbitals)
	printf("calculating n=%2d  l=%2d  :  Energy guess=%12.6lf\n", n, l, energy);
      
      atomLevels[nl].rho.resize(a.r_mesh.size());
      
      deepst_check_(&n, &l, &kappa, &energy, &a.vr(0,0), &a.r_mesh[0],
              &atomLevels[nl].rho[0], &a.h,
              &a.ztotss, &cLight, &nitmax, &tol, &a.jws, &last, &iter,
		    &last, &ipdeq, &found);
      // normalize the density:
      for(int ir=0; ir<a.r_mesh.size(); ir++)
        atomLevels[nl].rho[ir]=atomLevels[nl].rho[ir]; ///a.r_mesh[ir];
      integrateOneDim(a.r_mesh, atomLevels[nl].rho, rf);
      Real normalizationFactor = 1.0/rf[a.jws]; // rf[last-2];
      char fname[256];
      
      FILE *orbitalFile;
      if(writeOrbitals)
      {
	sprintf(fname,"orbital_density_n%d_l%d",n,l);
	orbitalFile=fopen(fname,"w");
      }
      
      for(int ir=0; ir<a.r_mesh.size(); ir++)
      {
        unnormalizedDensity[ir]=atomLevels[nl].rho[ir];
        atomLevels[nl].rho[ir]=atomLevels[nl].rho[ir]*normalizationFactor; // *a.r_mesh[ir];;
      }
      // test normalization
      integrateOneDim(a.r_mesh, atomLevels[nl].rho, rfNormalized);
      
      if(writeOrbitals)
      {
	for(int ir=0; ir<a.r_mesh.size(); ir++)
        {
	  fprintf(orbitalFile,"%4d %g  %g %g",ir,a.r_mesh[ir],unnormalizedDensity[ir],rf[ir]);
	  fprintf(orbitalFile," %g",atomLevels[nl].rho[ir]);
	  fprintf(orbitalFile," %g\n",rfNormalized[ir]);
	}
	fclose(orbitalFile);
      }
      
      atomLevels[nl].n = n;
      atomLevels[nl].l = l;
      atomLevels[nl].energy = energy;
      if(writeOrbitals)
	printf("n=%2d  l=%2d  :  Energy=%12.6lf\n", n, l, energy);
      nl++;
    }
  }
    // sort
  std::sort(atomLevels.begin(), atomLevels.end(), compareInitialAtomLevels());
//            [](myclass const & a, myclass const &b){return a.energy < b.energy;});
// fill core states:
  int coreTarget = a.zcorss + a.zsemss;
  int coreElectrons = 0;
  a.numc = 0;
  for (int i=0; i<atomLevels.size(); i++)
  {
    if (coreElectrons < coreTarget)
    {
      if(coreElectrons<a.zcorss)
      {
        printf("C ");
        atomLevels[i].type='C';
      } else {
        printf("S ");
        atomLevels[i].type='S';
      }
      coreElectrons += 2*(2*atomLevels[i].l+1);
      if(atomLevels[i].l == 0)
        a.numc += 1;
      else
        a.numc += 2;
    }
    else
    {
      printf("V ");
      atomLevels[i].type='V';
    }
    printf("n=%2d  l=%2d  :  Energy=%12.6lf\n",atomLevels[i].n,atomLevels[i].l,atomLevels[i].energy);
  }

  a.resizeCore(a.numc);
  int j = 0;
  for (int i=0; i<a.numc; i++)
  {
    a.ec(i,0) = a.ec(i,1) = atomLevels[j].energy;
    a.nc(i,0) = a.nc(i,1) = atomLevels[j].n;
    a.lc(i,0) = a.lc(i,1) = atomLevels[j].l;
    a.kc(i,0) = a.kc(i,1) = -atomLevels[j].l-1;
    if(atomLevels[j].type=='C') // accumulate the density for core levels
    {
      for(int ir=0; ir<a.r_mesh.size(); ir++)
      {
        a.corden(ir,0) += Real((3-a.nspin)*(-a.kc(i,0)))*atomLevels[j].rho[ir];
        a.corden(ir,1) += Real((3-a.nspin)*(-a.kc(i,1)))*atomLevels[j].rho[ir];
      }
    } else if(atomLevels[j].type=='S') { // accumulate the density for semi-core levels
      for(int ir=0; ir<a.r_mesh.size(); ir++)
      {
        a.semcor(ir,0) += Real((3-a.nspin)*(-a.kc(i,0)))*atomLevels[j].rho[ir];
        a.semcor(ir,1) += Real((3-a.nspin)*(-a.kc(i,1)))*atomLevels[j].rho[ir];
      }
    }
    if(atomLevels[j].l != 0)
    {
      i++;
      a.ec(i,0) = a.ec(i,1) = atomLevels[j].energy;
      a.nc(i,0) = a.nc(i,1) = atomLevels[j].n;
      a.lc(i,0) = a.lc(i,1) = atomLevels[j].l;
      a.kc(i,0) = a.kc(i,1) = atomLevels[j].l;
      if(atomLevels[j].type=='C') // accumulate the density for core levels
      {
        for(int ir=0; ir<a.r_mesh.size(); ir++)
        {
          a.corden(ir,0)+=Real((3-a.nspin)*a.kc(i,0))*atomLevels[j].rho[ir];
          a.corden(ir,1)+=Real((3-a.nspin)*a.kc(i,1))*atomLevels[j].rho[ir];
        }
      } else if(atomLevels[j].type=='S') { // accumulate the density for semi-core levels
        for(int ir=0; ir<a.r_mesh.size(); ir++)
        {
          a.semcor(ir,0)+=Real((3-a.nspin)*a.kc(i,0))*atomLevels[j].rho[ir];
          a.semcor(ir,1)+=Real((3-a.nspin)*a.kc(i,1))*atomLevels[j].rho[ir];
        }
      }
    }
    j++;
    // printf("%d %d\n",i,j);
  }
  // accumulate valence density
  for(int ir=0; ir<a.r_mesh.size(); ir++)
  {
    rhoValence(ir,0)=rhoValence(ir,1)=0.0;
  }
  int numValenceRemaining = a.zvalss;
  j=0; while(atomLevels[j].type != 'V') j++;
  while(numValenceRemaining>0)
  {
    int multiplicity = 2*atomLevels[j].l+1;
    if(2*multiplicity <= numValenceRemaining)
    {
      for(int ir=0; ir<a.r_mesh.size(); ir++)
      {
	rhoValence(ir,0) += Real(multiplicity*(3-a.nspin))*atomLevels[j].rho[ir];
	rhoValence(ir,1) += Real(multiplicity*(3-a.nspin))*atomLevels[j].rho[ir];
      }
      numValenceRemaining -= 2*multiplicity;
      j++;
    } else {
      for(int ir=0; ir<a.r_mesh.size(); ir++)
      {
	rhoValence(ir,0) += 0.5*Real(numValenceRemaining*(3-a.nspin))*atomLevels[j].rho[ir];
	rhoValence(ir,1) += 0.5*Real(numValenceRemaining*(3-a.nspin))*atomLevels[j].rho[ir];
      }
      numValenceRemaining = 0;
      j++;
    }
  }
  
  return coreElectrons;
}

void initializeAtom(AtomData &a)
{
  a.generateRadialMesh();

  // for analysis purposes gnerate gnuplot files for the atom
  bool generatePlot=true;
  FILE * plotFile;
  if(generatePlot) plotFile=fopen("initializeAtom.plot","w");
  
  // inititalize potential to be V(r) = -2Z/r
  // add a potential corresponding to a gaussian charge distribution
  // for the core charge with \sigma=0.2*rmt 
  // (vr stores V(r) * r)
  Real sigmaSqrt2Inv = 1.0 / (0.2 * std::sqrt(2.0) * a.rmt);
  Real q = a.zcorss + a.zsemss;
  for (int ir=0; ir<a.r_mesh.size(); ir++)
  {
    a.vr(ir,0) = a.vr(ir,1) = -2.0*a.ztotss; // +q*approxErfc(a.r_mesh[ir]*sigmaSqrt2Inv);
    // clear core and semi-core densities
    a.corden(ir,0) = a.corden(ir,1) = a.semcor(ir,0) = a.semcor(ir,1) = 0.0;
    a.rhotot(ir,0) = a.rhotot(ir,1) = a.rhoNew(ir,0) = a.rhoNew(ir,1) = 0.0;
  }
  // add homogeneous charge density inside the sphere from valence electrons
  // do make the site neutral. (i.e. vr(r_mesh.size(),*)=0)
  // q=-a.vr(a.r_mesh.size()-1,0)/a.r_mesh[a.r_mesh.size()-1];
  // for(int  ir=0; ir<a.r_mesh.size(); ir++)
  // {
  //   Real r=a.r_mesh[ir];
  //   a.vr(ir,0)+=q*r*r;
  //   a.vr(ir,1)+=q*r*r;
  // }

  // find the lowest zcorss+zsemss levels
  // maximum number of states to test: (lmax+1)*(lmax+2)/2
  // assuming kappa degeneracy, iterating over principal quantum number and l
  // we treat core and semi-core states the same for the initialization,
  // using deepst for both

  int lmaxCore = std::min(a.lmax, 3);   // only consider s,p,d,f electrons (we should never need g or higher)

  int  numAtomLevels = ((lmaxCore+1) * (lmaxCore+2)) / 2;
  int  kappa         = 0;
  Real energy        = 0.0;
  Real cLight        = 2.0 * 137.0359895;
  int  nitmax        = 100;                // maximum number of iterations
  Real tol           = 1.0e-8;             // energy tolerance
  int  last          = a.r_mesh.size();
  int  ipdeq         = 5;
  int  iter          = 0;
  std::vector<InitialAtomLevels> atomLevels(numAtomLevels);
  Matrix<Real> rhoValence(a.r_mesh.size(),2);
  
  int coreTarget = a.zcorss + a.zsemss;
  int coreElectrons = calculateAtomLevels(a, atomLevels, rhoValence, true);
  
  if (coreElectrons != coreTarget)
  {
    printf("Warning: initializeAtom can't satisfy the core electron requirement:\n  Target: %d (%lf + %lf)\n  Actual: %d\n",
           coreTarget, a.zcorss, a.zsemss, coreElectrons);
  }

  
  // add charge density contributions
  // take constant charge density from valence electrons
  // (add small spliting in spin up/down valence densities to initialize magnetic calculations)
  Real delta=0.1;
  
  // Real valenceDensityUp, valenceDensityDown;
  Real vmt, vmt1;
  vmt=0.0;
  vmt1=0.0;
  Real atomVolume=(4.0/3.0)*M_PI*a.rws*a.rws*a.rws;  // a.rmt*a.rmt*a.rmt;
  /*
  if(a.nspin==1)
  {
    for(int ir=0; ir<a.r_mesh.size(); ir++)
      rhoValence(ir,0)=rhoValence(ir,1)=a.zvalss/atomVolume;
  } else {
    for(int ir=0; ir<a.r_mesh.size(); ir++)
    {
      rhoValence(ir,0)=(0.5*a.zvalss+delta)/atomVolume;
      rhoValence(ir,1)=(0.5*a.zvalss-delta)/atomVolume;
    }
  }
  */
  for(int ir=0; ir<a.r_mesh.size(); ir++)
  {
    a.rhotot(ir,0)=a.rhoNew(ir,0)=a.corden(ir,0)+a.semcor(ir,0)+rhoValence(ir,0);
    a.rhotot(ir,1)=a.rhoNew(ir,1)=a.corden(ir,1)+a.semcor(ir,1)+rhoValence(ir,1);
  }
  // calculate potential from initial density guess
  std::vector<Real> chargeDensity(a.r_mesh.size()+1);
  std::vector<Real> electrostaticPotential(a.r_mesh.size()+1);
  // std::vector<Real> exchangeCorrelationPotential(a.r_mesh.size()+1);
  std::vector<Real> mesh0(a.r_mesh.size()+1);

  // printf("jmt =%d\n",a.jmt);

  Matrix<Real> rhoTemp(a.jmt+2,2); 
  std::vector<Real> rTemp(a.jmt+3);
  rTemp[0] = 0.0;
  for (int j=0; j<a.jmt+2; j++)
  {
    //rTemp's indices need to be shifted by 1 for passing into getqm_mt!
    rTemp[j+1] = std::sqrt(a.r_mesh[j]);
  }

  
  chargeDensity[0]=0.0;
  mesh0[0]=0.0;
  if(a.nspin==1)
    for(int ir=0; ir<a.r_mesh.size(); ir++)
    {
      chargeDensity[ir+1]=a.rhotot(ir,0); // factor 2 for solution of radial poisson eq. : v(r) = 8 pi/r \int_0^r r'^2 \rho(r')dr' 
      mesh0[ir+1]=a.r_mesh[ir];
    }
  else
    for(int ir=0; ir<a.r_mesh.size(); ir++)
    {
      chargeDensity[ir+1]=2.0*(a.rhotot(ir,0)+a.rhotot(ir,1));
      mesh0[ir+1]=a.r_mesh[ir];
    }
  
  integrateOneDim(mesh0, chargeDensity, electrostaticPotential);

  // calculate the exchange correlation potential
  for(int is=0; is<a.nspin; is++)
  {
    a.qInt = 0.0; // a.qtotws - a.qtotmt;
    Real spin = 1.0 - Real(is) * 2.0;
    Real ro3 = std::pow( (4.0*M_PI/3.0) * a.rhoInt, -1.0/3.0);
    Real dz = 0.0; //a.mInt / a.qInt;
    int iexch = 0; // von Barth-Hedin
    int mtasa = 0;
    int jmt = a.jmt;
    
    newexchg_(&a.nspin, &spin, &a.rhotot(0,0), &a.rhotot(0,a.nspin-1),
	      &a.exchangeCorrelationPotential(0,is), &a.exchangeCorrelationEnergy(0,is),
	      &a.exchangeCorrelationV[is], &a.exchangeCorrelationE,
	      &ro3, &dz, &a.r_mesh[0], &a.jmt, &iexch);

    newpot_(&a.nspin, &a.ztotss, &a.rhotot(0,0),
	    &a.rhotot(0,a.nspin-1), &rhoTemp(0,0),
	    &a.vr(0,is), &a.vrNew(0,is), &a.vrms[is],
	    &a.exchangeCorrelationPotential(0,is), &vmt1, &vmt, &a.exchangeCorrelationV[is], &rTemp[0],
	    &a.jmt, &a.rInscribed, &a.rInscribed,
	    &mtasa, &iexch);
  }
  
  Real potMix = 0.1; // 0.1
  for(int ir=0; ir<a.r_mesh.size(); ir++)
  {
    // a.vr(ir,0)= potMix*(-2.0*a.ztotss+electrostaticPotential[ir+1]+a.exchangeCorrelationPotential(ir,0))+(1.0-potMix)*a.vr(ir,0);
    // a.vr(ir,1)= potMix*(-2.0*a.ztotss+electrostaticPotential[ir+1]+a.exchangeCorrelationPotential(ir,1))+(1.0-potMix)*a.vr(ir,1);
    a.vr(ir,0)= potMix*a.vrNew(ir,0)+(1.0-potMix)*a.vr(ir,0);
    a.vr(ir,1)= potMix*a.vrNew(ir,1)+(1.0-potMix)*a.vr(ir,1);
  }
  for(int ir=a.jmt-1; ir<a.r_mesh.size(); ir++)
  {
    a.vr(ir,0)= a.vr(a.jmt-2,0);
    a.vr(ir,1)= a.vr(a.jmt-2,1);
  }
  if(generatePlot)
  {
    fprintf(plotFile,"set term pdf\nset outp 'initialPotential_start.pdf'\n");
    fprintf(plotFile,"set xrange [0:%lf]\n",a.r_mesh[a.r_mesh.size()-1]);

    fprintf(plotFile,"plot '-' with lines title 'Vr spin up'\n");
    for(int ir=0; ir<a.r_mesh.size(); ir++)
    {
      fprintf(plotFile,"%18.12lf   %18.12lf\n",a.r_mesh[ir],a.vr(ir,0));
    }
    fprintf(plotFile,"e\n");
    
    fprintf(plotFile,"plot '-' with lines title 'Vr spin down'\n");
    for(int ir=0; ir<a.r_mesh.size(); ir++)
    {
      fprintf(plotFile,"%18.12lf   %18.12lf\n",a.r_mesh[ir],a.vr(ir,1));
    }
    fprintf(plotFile,"e\n");
    
    fprintf(plotFile,"plot '-' with lines title 'rho spin up'\n");
    for(int ir=0; ir<a.r_mesh.size(); ir++)
    {
      fprintf(plotFile,"%18.12lf   %18.12lf\n",a.r_mesh[ir],a.rhotot(ir,0));
    }
    fprintf(plotFile,"e\n");
    
    fprintf(plotFile,"plot '-' with lines title 'rho spin down'\n");
    for(int ir=0; ir<a.r_mesh.size(); ir++)
    {
      fprintf(plotFile,"%18.12lf   %18.12lf\n",a.r_mesh[ir],a.rhotot(ir,1));
    }
    fprintf(plotFile,"e\n");
  }
  
  potMix = 0.02;
  
// iteration for preliminary potential
  for(int iteration=0; iteration<20; iteration++)
  {
    // printf("iter %d\n",iteration);
    for(int ir=0; ir<a.r_mesh.size(); ir++)
    {
      a.rhotot(ir,0)=a.rhoNew(ir,0)=0.0;
      a.rhotot(ir,1)=a.rhoNew(ir,1)=0.0;
    }

    calculateAtomLevels(a, atomLevels, rhoValence);

    for(int ir=0; ir<a.r_mesh.size(); ir++)
    {
      a.rhotot(ir,0)=a.rhoNew(ir,0)=a.corden(ir,0)+a.semcor(ir,0)+rhoValence(ir,0);
      a.rhotot(ir,1)=a.rhoNew(ir,1)=a.corden(ir,1)+a.semcor(ir,1)+rhoValence(ir,1);
    }

    chargeDensity[0]=0.0;
    mesh0[0]=0.0;
    if(a.nspin==1)
      for(int ir=0; ir<a.r_mesh.size(); ir++)
      {
        chargeDensity[ir+1]=2.0*a.rhotot(ir,0); // factor 2 for solution of radial poisson eq. : v(r) = 8 pi/r \int_0^r r'^2 \rho(r')dr' 
        mesh0[ir+1]=a.r_mesh[ir];
      }
    else
      for(int ir=0; ir<a.r_mesh.size(); ir++)
      {
        chargeDensity[ir+1]=2.0*(a.rhotot(ir,0)+a.rhotot(ir,1));
        mesh0[ir+1]=a.r_mesh[ir];
      }
    integrateOneDim(mesh0, chargeDensity, electrostaticPotential);

  // calculate the exchange correlation potential
    for(int is=0; is<a.nspin; is++)
    {
      a.qInt = 0.0; // a.qtotws - a.qtotmt;
      Real spin = 1.0 - Real(is) * 2.0;
      Real ro3 = std::pow( (4.0*M_PI/3.0) * a.rhoInt, -1.0/3.0);
      Real dz = 0.0; //a.mInt / a.qInt;
      int iexch = 0; // von Barth-Hedin
      int mtasa=0;
      int jmt = a.jmt;
      
      newexchg_(&a.nspin, &spin, &a.rhotot(0,0), &a.rhotot(0,a.nspin-1),
		&a.exchangeCorrelationPotential(0,is), &a.exchangeCorrelationEnergy(0,is),
		&a.exchangeCorrelationV[is], &a.exchangeCorrelationE,
		&ro3, &dz, &a.r_mesh[0], &a.jmt, &iexch);

      
      newpot_(&a.nspin, &a.ztotss, &a.rhotot(0,0),
	      &a.rhotot(0,a.nspin-1), &rhoTemp(0,0),
	      &a.vr(0,is), &a.vrNew(0,is), &a.vrms[is],
	      &a.exchangeCorrelationPotential(0,is), &vmt1, &vmt, &a.exchangeCorrelationV[is], &rTemp[0],
	      &a.jmt, &a.rInscribed, &a.rInscribed,
	      &mtasa, &iexch);
    }
    
    potMix = 0.1; // 0.1
    for(int ir=0; ir<a.r_mesh.size(); ir++)
    {
      // a.vr(ir,0)= potMix*(-2.0*a.ztotss+electrostaticPotential[ir+1]+a.exchangeCorrelationPotential(ir,0))+(1.0-potMix)*a.vr(ir,0);
      // a.vr(ir,1)= potMix*(-2.0*a.ztotss+electrostaticPotential[ir+1]+a.exchangeCorrelationPotential(ir,1))+(1.0-potMix)*a.vr(ir,1);
      a.vr(ir,0)= potMix*a.vrNew(ir,0)+(1.0-potMix)*a.vr(ir,0);
      a.vr(ir,1)= potMix*a.vrNew(ir,1)+(1.0-potMix)*a.vr(ir,1);
    }
    for(int ir=a.jmt-1; ir<a.r_mesh.size(); ir++)
    {
      a.vr(ir,0)= a.vr(a.jmt-2,0);
      a.vr(ir,1)= a.vr(a.jmt-2,1);
    }
  }

  if(generatePlot)
  {
    fprintf(plotFile,"set term pdf\nset outp 'initialPotential_final.pdf'\n");
    fprintf(plotFile,"set xrange [0:%lf]\n",a.r_mesh[a.r_mesh.size()-1]);

    fprintf(plotFile,"plot '-' with lines title 'Vr spin up'\n");
    for(int ir=0; ir<a.r_mesh.size(); ir++)
    {
      fprintf(plotFile,"%18.12lf   %18.12lf\n",a.r_mesh[ir],a.vr(ir,0));
    }
    fprintf(plotFile,"e\n");
    
    fprintf(plotFile,"plot '-' with lines title 'Vr spin down'\n");
    for(int ir=0; ir<a.r_mesh.size(); ir++)
    {
      fprintf(plotFile,"%18.12lf   %18.12lf\n",a.r_mesh[ir],a.vr(ir,1));
    }
    fprintf(plotFile,"e\n");
    
    fprintf(plotFile,"plot '-' with lines title 'rho spin up'\n");
    for(int ir=0; ir<a.r_mesh.size(); ir++)
    {
      fprintf(plotFile,"%18.12lf   %18.12lf\n",a.r_mesh[ir],a.rhotot(ir,0));
    }
    fprintf(plotFile,"e\n");
    
    fprintf(plotFile,"plot '-' with lines title 'rho spin down'\n");
    for(int ir=0; ir<a.r_mesh.size(); ir++)
    {
      fprintf(plotFile,"%18.12lf   %18.12lf\n",a.r_mesh[ir],a.rhotot(ir,1));
    }
    fprintf(plotFile,"e\n");
  }
  
  if(generatePlot) fclose(plotFile);
  // exit(1);
}


int initializeNewPotentials(LSMSCommunication &comm,LSMSSystemParameters &lsms, CrystalParameters &crystal, LocalTypeInfo &local)
{
  for (int i=0; i<local.num_local; i++)
  {
    printf("Initializing potential %d.\n", i);
    snprintf(local.atom[i].header, 40, "Potential Initialized by LSMS_3");
    for (int j=31; j<80; j++)
      local.atom[i].header[j] = ' ';
    local.atom[i].resizePotential(lsms.global.iprpts);
    local.atom[i].ztotss = (Real)crystal.types[local.global_id[i]].Z;
    local.atom[i].zcorss = (Real)crystal.types[local.global_id[i]].Zc;
    local.atom[i].zsemss = (Real)crystal.types[local.global_id[i]].Zs;
    local.atom[i].zvalss = (Real)crystal.types[local.global_id[i]].Zv;
    local.atom[i].vdif = 0.0;
    local.atom[i].vdifNew = 0.0;
    local.atom[i].spinFlipped = false;
    local.atom[i].xvalws[0] = local.atom[i].xvalws[1] = 0.5 * local.atom[i].zvalss;
    local.atom[i].lmax = crystal.types[local.global_id[i]].lmax;
    local.atom[i].kkrsz = (local.atom[i].lmax + 1) * (local.atom[i].lmax + 1);

    local.atom[i].xstart = -11.1309674;
    local.atom[i].jmt = 1001;
    local.atom[i].jws = 1012;
    printf("rmt=%lf\n",local.atom[i].rmt);
    local.atom[i].generateRadialMesh();
    local.atom[i].nspin = 2;
    local.atom[i].evec[0] = local.atom[i].evec[1] = 0.0;
    local.atom[i].evec[2] = 1.0;

    initializeAtom(local.atom[i]);
    local.atom[i].alat = local.atom[i].rmt;
    local.atom[i].efermi = 0.5;
  }

  lsms.chempot = 0.5;
  calculateCoreStates(comm, lsms, local);
  return 0;
}

int initializeNewAlloyBank(LSMSCommunication &comm, LSMSSystemParameters &lsms, AlloyMixingDesc &alloyDesc, AlloyAtomBank &alloyBank)
{
  
  printf("Initializing new alloy bank\n");

  for(int id = 0, i = 0; i < alloyBank.size(); i++)
  for(int j = 0; j < alloyBank[i].size(); j++, id++) {

    AtomData &atom = alloyBank[i][j];

    printf("Initializing potential %d.\n",id);
    snprintf(atom.header, 40, "Potential Initialized by LSMS_3");

    for(int j=31;j<80;j++) atom.header[j]=' ';
    atom.resizePotential(lsms.global.iprpts);
    atom.ztotss=(Real)alloyDesc[i][j].Z;
    atom.zcorss=(Real)alloyDesc[i][j].Zc;
    atom.zsemss=(Real)alloyDesc[i][j].Zs;
    atom.zvalss=(Real)alloyDesc[i][j].Zv;
    atom.alloy_class=(Real)alloyDesc[i][j].alloy_class;
    atom.vdif=0.0;
    atom.xvalws[0]=atom.xvalws[1]=0.5*atom.zvalss;
    atom.lmax=alloyDesc[i][j].lmax;
    atom.kkrsz=(atom.lmax+1)*(atom.lmax+1);

    atom.xstart=-11.1309674;
    atom.jmt=1001;
    atom.jws=1001;
    atom.nspin=2;
    atom.evec[0]=atom.evec[1]=0.0; atom.evec[2]=1.0;

    initializeAtom(atom);
    atom.alat=atom.rmt;
    atom.efermi=0.5;
  }
  calculateCoreStates(comm, lsms, alloyBank);
  return 0;
}
