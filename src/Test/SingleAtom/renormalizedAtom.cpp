// Renormalized Atom

#include <tuple>
#include <algorithm>

// typedef double Real;
#include "Real.hpp"
#include "Complex.hpp"
#include "Matrix.hpp"

#include "../../Misc/bulirschStoerIntegrator.hpp"
#include "../../Misc/rationalFit.hpp"
#include "../../Misc/integrateOneDim.cpp"

/*
// integrate a system of n coupled ordinary differential equations
// from x0 to x1 with initial values in y0
// return values at x1 in y1
// rhs is a function that can evaluate the right hand side of the system of ODEs
// anywhere in the interval [x0,x1]
// eps is the desired max error
// returns an error code: 0 -> success
// templated on the type of real numbers
template <typename Rx, typename Ry>
int bulirschStoerIntegrator(Rx x0, Rx x1, Ry *y0,
                            Ry *y1, int n, std::function<void(Rx x, Ry* y, Ry* dy)> rhs,
                            Rx eps=1.0e-12)

/// Given a table of function values r[i] -> f(r[i])
/// find the interpolated value of f(x)
template<typename T>
T interpolate(std::vector<T> &r, std::vector<T> &f, T x)
*/

// Dirac RHS for radial electrostatic potential without magnetic field:
// in: r, pq[0] = p (= r*g) , pd[1] = q = (= r*f)
// out: dpq[0] = dp/dr , dpq[1] = dq/dr

void diracRHS(Real r, Real *pq, Real *dpq, std::vector<Real> &r_mesh, std::vector<Real> &vr,
              Real e, Real kappa)
{
  const Real cLight = 274.072;
  Real v = interpolate<Real>(r_mesh, vr, r)/r;
  // note: e = W - mc^2
  // dP/dr = -kappa/r * P(r) + 1/(c hbar) * (e - V(r) + 2mc^2) * Q(r)
  // dQ/dr =  kappa/r * Q(r) - 1/(c hbar) * (e - V(r)        ) * P(r)
  dpq[0] = -(kappa/r) * pq[0] + ((e-v)/cLight + cLight) * pq[1];
  dpq[1] =  (kappa/r) * pq[1] - ((e-v)/cLight         ) * pq[0];
}

template<typename PQType>
void diracBoundaryConditionNearOrigin(Real r, PQType *pq, Real Z, Real kappa)
{
  const Real cLight = 274.072;
  Real zeta = 2.0*Z/cLight; // Z * alpha

  pq[0] = 1.0e-20;
  // pq[1] = ((kappa + std::sqrt(kappa*kappa - zeta*zeta))/zeta)*cLight * pq[0];
  pq[1] = ((kappa + std::sqrt(kappa*kappa - zeta*zeta))/zeta) * pq[0];
}

Real energyGuess(int n, Real Z, Real kappa)
{
  const Real cLight = 274.072;
  Real zeta = 2.0*Z/cLight; // Z * alpha

//  return (2.0*Z*Z)*(-1.0/(2.0*Real(n*n)));

  return (2.0*Z*Z)*(-1.0/(2.0*Real(n*n)) +
                      (zeta*zeta)*(3.0/(8.0*Real(n*n*n*n)) - 1.0/(2.0*Real(n*n*n)*std::abs(kappa)))); 
}

void generateRadialMesh(std::vector<Real> &r_mesh, int N, Real r0, Real rN)
{
    if (N != r_mesh.size()) r_mesh.resize(N);
    Real x0 = std::log(r0);
    Real xN = std::log(rN);
    Real h = (xN-x0) / (N-1);
    for(int j=0; j<N; j++)
    {
      Real x_mesh = x0 + (Real)j*h;
      r_mesh[j] = std::exp(x_mesh);
    }
}

int countNodes(Matrix<Real> &y, int c=0)
{
  int n=0;
  bool s = std::signbit(y(c,0));
  for(int i=1; i<y.n_col()-1; i++)
  {
    if(s != std::signbit(y(c,i)))
    {
      n++;
      s = std::signbit(y(c,i));
    }
  }
  return n;
}

void integrateDirac(std::vector<Real> &r_mesh, Real atomicNumber, std::vector<Real> &vr,
                    Real energy, Real kappa, Matrix<Real> &pq)
{
  // printf("integrateDirac: energy=%lg\n",energy);
// outward integration
  diracBoundaryConditionNearOrigin(r_mesh[0], &pq(0,0), Real(atomicNumber), kappa);

  // for(int i=1; i<iFitting; i++)
  for(int i=1; i<r_mesh.size(); i++)
  {
    if(bulirschStoerIntegrator<Real, Real>(r_mesh[i-1], r_mesh[i], &pq(0,i-1), &pq(0,i), 2,
                           [&](Real r, Real *y, Real *dy)
                           {diracRHS(r, y, dy, r_mesh, vr, energy, kappa);},
                           1.0e-12))
    {
      printf("integration did not succeed: %d:%lg -> %d:%lg!\n",i-1,r_mesh[i-1],i,r_mesh[i]);
    }
  }
}

void calculateRadialDensity(std::vector<Real> &r_mesh, Matrix<Real> &pq, std::vector<Real> &rho)
{
  std::vector<Real> rhoIntegrated(r_mesh.size());

  for(int i=0; i<r_mesh.size(); i++)
  {
    rho[i] = std::abs(pq(0,i))*std::abs(pq(0,i))
      + std::abs(pq(1,i))*std::abs(pq(1,i));
  }

  // integrateOneDimSpherical(r_mesh, rho, rhoIntegrated);
  // the p & q already contain a factor of r,
  // thus the rho above is actually r^2 rho and a standard integration
  // is to be performed to yield the spherical integral of |psi|^2
  // we only need a facto of 4 pi
  integrateOneDim(r_mesh, rho, rhoIntegrated);

  Real normalization = 1.0/(4.0 * M_PI * rhoIntegrated[rhoIntegrated.size()-1]);

  for(int i=0; i<r_mesh.size(); i++)
  {
    rho[i] *= normalization;
  }
}

Real findDiracEigenvalue(std::vector<Real> &r_mesh, std::vector<Real> &vr, Real atomicNumber,
			 int principalQuantumNumber, int kappa, Real energy,
			 Matrix<Real> &pq)
{
  Real energyUpper, energyLower;
  int l = kappa;
  if(l<0) l = -kappa - 1;

  int targetNumberOfNodes = principalQuantumNumber - l - 1;

// find energy bounds

  integrateDirac(r_mesh, atomicNumber, vr, energy, kappa, pq);

  int numNodesP = countNodes(pq,0);
  int numNodesQ = countNodes(pq,1);

  if(numNodesP>targetNumberOfNodes)
  {
    energyUpper = energy;
    while(numNodesP>targetNumberOfNodes)
    {
      energy -= 1.0;
      integrateDirac(r_mesh, atomicNumber, vr, energy, kappa, pq);
      numNodesP = countNodes(pq,0);
    }
    energyLower = energy;
  } else {
    energyLower = energy;
    while(numNodesP<=targetNumberOfNodes)
    {
      energy += 1.0;
      integrateDirac(r_mesh, atomicNumber, vr, energy, kappa, pq);
      numNodesP = countNodes(pq,0);
    }
    energyUpper = energy;
  }

  Real energyEpsilon = 1.0e-15;
  while(std::abs((energyUpper - energyLower)/energy) > energyEpsilon)
  {
    energy = energyLower + 0.5*(energyUpper-energyLower);
    integrateDirac(r_mesh, atomicNumber, vr, energy, kappa, pq);
    numNodesP = countNodes(pq,0);
    if(numNodesP>targetNumberOfNodes)
    {
      energyUpper = energy;
    } else {
      energyLower = energy;
    }
    // printf("%lf %lf  %lg\n",energyLower, energyUpper, energyUpper - energyLower);
  }

  return energy;
}

class AtomOrbital {
public:
  int n,kappa;
  Real energy;
  std::vector<Real> rho; // radial charge density contribution from an electron in this level
  // char type; // 'C' for core lectron, 'S' for semi-core and 'V' for valence electrons
} ;


// generate a list of (n, kappa) orbitals
void initOrbitals(int Z, std::vector<std::tuple<int, int>> &orbitals)
{
  orbitals.clear();
  // He: 1s^2
  orbitals.push_back(std::make_tuple<int,int>(1,-1));
  if(Z<3) return;
  // Ne: [He] 2s^2 2p^6
  orbitals.push_back(std::make_tuple<int,int>(2,-1));
  orbitals.push_back(std::make_tuple<int,int>(2, 1));
  orbitals.push_back(std::make_tuple<int,int>(2,-2));
  if(Z<11) return;
  // Ar: [Ne] 3s^2 3p^6
  orbitals.push_back(std::make_tuple<int,int>(3,-1));
  orbitals.push_back(std::make_tuple<int,int>(3, 1));
  orbitals.push_back(std::make_tuple<int,int>(3,-2));
  if(Z<19) return;
  // Kr: [Ar] 3d^10 4s^2 4p^6
  orbitals.push_back(std::make_tuple<int,int>(3, 2));
  orbitals.push_back(std::make_tuple<int,int>(3,-3));
  orbitals.push_back(std::make_tuple<int,int>(4,-1));
  orbitals.push_back(std::make_tuple<int,int>(4, 1));
  orbitals.push_back(std::make_tuple<int,int>(4,-2));
  if(Z<37) return;
  // Xe: [Kr] 4d^10 5s^2 5p^6
  orbitals.push_back(std::make_tuple<int,int>(4, 2));
  orbitals.push_back(std::make_tuple<int,int>(4,-3));
  orbitals.push_back(std::make_tuple<int,int>(5,-1));
  orbitals.push_back(std::make_tuple<int,int>(5, 1));
  orbitals.push_back(std::make_tuple<int,int>(5,-2));
  if(Z<55) return;
  // Rn: [Xe] 4f^14 5d^10 6s^2 6p^6
  orbitals.push_back(std::make_tuple<int,int>(4, 3));
  orbitals.push_back(std::make_tuple<int,int>(4,-4));
  orbitals.push_back(std::make_tuple<int,int>(5, 2));
  orbitals.push_back(std::make_tuple<int,int>(5,-3));
  orbitals.push_back(std::make_tuple<int,int>(6,-1));
  orbitals.push_back(std::make_tuple<int,int>(6, 1));
  orbitals.push_back(std::make_tuple<int,int>(6,-2));
  if(Z<87) return;
  // Og: [Rn] 5f^14 6d^10 7s^2 7p^6
  orbitals.push_back(std::make_tuple<int,int>(5, 3));
  orbitals.push_back(std::make_tuple<int,int>(5,-4));
  orbitals.push_back(std::make_tuple<int,int>(6, 2));
  orbitals.push_back(std::make_tuple<int,int>(6,-3));
  orbitals.push_back(std::make_tuple<int,int>(7,-1));
  orbitals.push_back(std::make_tuple<int,int>(7, 1));
  orbitals.push_back(std::make_tuple<int,int>(7,-2));
}

void sphericalPoisson(Real *rho, Real *r, int nr, Real *vr)
{
  
}

void printUsage(const char *name, Real R)
{
  printf("Usage: %s Z [R]\n",name);
  printf("       Z: atomic number\n");
  printf("       R: atomic sphere radius (optional, default=%lf)\n",R);
}

int main(int argc, char *argv[])
{
  int atomicNumber = 29; // test copper
  int principalQuantumNumber = 1;
  int kappa = -1; // l=0
  int atomRadius = 3.0;

  if(argc != 2 && argc != 3)
  {
    printUsage(argv[0], atomRadius);
    exit(1);
  }
  atomicNumber = atoi(argv[1]);
  if(argc == 3)
    atomRadius = atof(argv[2]);

  int maxPrincipalQuantumNumber;

  std::vector<std::tuple<int, int>> orbitals;
  initOrbitals(atomicNumber, orbitals);

  printf("number of orbitals to be computed: %d\n",orbitals.size());
  
  std::vector<Real> r_mesh, vr;
  std::vector<Matrix<Real> > pqs;
  std::vector<AtomOrbital> orbitalEnergiesAndDensities;
  pqs.resize(orbitals.size());
  orbitalEnergiesAndDensities.resize(orbitals.size());

// initialize r_mesh (unit of length is the Bohr radius)
  generateRadialMesh(r_mesh, 1500, 1.5e-10, atomRadius);

// initialize Z/r potential e^2=2
  vr.resize(r_mesh.size());
  for(int i=0; i<r_mesh.size(); i++) vr[i] = -2.0*Real(atomicNumber);

  printf(" n  kappa  energy\n");
  for(int i=0; i<orbitals.size(); i++)
  {
    Matrix<Real> &pq = pqs[i];
    principalQuantumNumber = std::get<0>(orbitals[i]);
    kappa = std::get<1>(orbitals[i]);
    pq.resize(2,r_mesh.size());
  
    Real energy = energyGuess(principalQuantumNumber, atomicNumber, kappa);
  // printf("# energy: %lg Ry\n",energy);

    orbitalEnergiesAndDensities[i].n = principalQuantumNumber;
    orbitalEnergiesAndDensities[i].kappa = kappa;
    orbitalEnergiesAndDensities[i].energy = findDiracEigenvalue(r_mesh, vr, atomicNumber,
					     principalQuantumNumber, kappa, energy,
					     pq);

    printf(" %d    %2d    %lg Ry\n",principalQuantumNumber, kappa,orbitalEnergiesAndDensities[i].energy);

    orbitalEnergiesAndDensities[i].rho.resize(r_mesh.size());
    calculateRadialDensity(r_mesh, pq, orbitalEnergiesAndDensities[i].rho);
    
    int l = kappa;
    if(l<0) l = -kappa - 1;
    int targetNumberOfNodes = principalQuantumNumber - l - 1;
    int numNodesP = countNodes(pq,0);
    int numNodesQ = countNodes(pq,1);
  }
  // accumulate the charge from all the electrons in the atom
  std::vector<Real> rhotot;
  rhotot.resize(r_mesh.size());
  int electronsMissing = atomicNumber; // still need Z electrons
  std::sort(orbitalEnergiesAndDensities.begin(), orbitalEnergiesAndDensities.end(),
	    [](AtomOrbital const & a, AtomOrbital const &b){return a.energy < b.energy;});
  int orbitalIdx=0; // the index of the current orbital
  while(electronsMissing>0)
  {
    if(electronsMissing >= 2*std::abs(orbitalEnergiesAndDensities[orbitalIdx].kappa)) // a filled orbital
    {
      for(int i=0; i<rhotot.size(); i++)
	rhotot[i] += orbitalEnergiesAndDensities[orbitalIdx].rho[i]
	  * 2 * std::abs(orbitalEnergiesAndDensities[orbitalIdx].kappa);
      electronsMissing -= 2 * std::abs(orbitalEnergiesAndDensities[orbitalIdx].kappa);
      printf("%d %2d %lg Ry: filled (%2d electrons)\n",
	     orbitalEnergiesAndDensities[orbitalIdx].n,
	     orbitalEnergiesAndDensities[orbitalIdx].kappa,
	     orbitalEnergiesAndDensities[orbitalIdx].energy,
	     2 * std::abs(orbitalEnergiesAndDensities[orbitalIdx].kappa));
    } else {
      for(int i=0; i<rhotot.size(); i++)
	rhotot[i] += orbitalEnergiesAndDensities[orbitalIdx].rho[i]
	  * electronsMissing;
      
      printf("%d %2d %lg Ry: partially filled (%2d electrons)\n",
	     orbitalEnergiesAndDensities[orbitalIdx].n,
	     orbitalEnergiesAndDensities[orbitalIdx].kappa,
	     orbitalEnergiesAndDensities[orbitalIdx].energy,
	     electronsMissing);

      electronsMissing = 0;
    }
    orbitalIdx++;
  }

  /*
  printf("# energy: %lg Ry\n",energy);
  printf("# nodes: P:%d Q:%d\n",numNodesP, numNodesQ);
  printf("# target number of nodes: %d\n",targetNumberOfNodes);

  for(int i=0; i<r_mesh.size(); i++)
  {
    Real rho = std::abs(pq(0,i))*std::abs(pq(0,i))
             + std::abs(pq(1,i))*std::abs(pq(1,i));
    printf("%5d %lg %lg %lg %lg\n",i, r_mesh[i], rho, pq(0,i), pq(1,i));
  }
  */

  return 0;
}
