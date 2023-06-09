/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */

#include "setupVorpol.hpp"

#include <cstdio>

void setupVorpol(LSMSSystemParameters &lsms, CrystalParameters &crystal, LocalTypeInfo &local)
{

  const bool testingVorpol = false;

  int ipvp=50; // parameter from LSMS_1.9 vorpol.h
  int iprcrit=260;
  int ipnode=ipvp*(ipvp-1);
  int ipcorn=(ipvp*(ipvp-1)*(ipvp-2))/6;
  int ipedge=(ipvp*ipvp-1)/2;
  const Real sphereVolumeFactor=4.0*M_PI/3.0;

  const int numSwitchAlgorithms = 500;
  int maxClusterSize = 0;

  double timeSetupVorpol = MPI_Wtime();

  Array3d<Real> vplanes(3,ipvp, local.num_local);
  std::vector<int> nvplanes(local.num_local);

  Array3d<Real> vplanesTest(3,ipvp, local.num_local);
  std::vector<int> nvplanesTest(local.num_local);

  bool useOldAlgorithm = false;

  if(lsms.global.iprint >= 0)
  {
    if(useOldAlgorithm)
     printf("Constructing Voronoi polyhedra using the OLD algorithm.\n");
    else
     printf("Constructing Voronoi polyhedra using the NEW algorithm.\n");
  }


  for(int i=0; i<local.num_local; i++)
  {
    if(local.atom[i].vpClusterGlobalIdx.size() > maxClusterSize)
      maxClusterSize = local.atom[i].vpClusterGlobalIdx.size();

    if(local.atom[i].vpClusterGlobalIdx.size() == 0)
      useOldAlgorithm = true;
  }
//  if(crystal.num_atoms < numSwitchAlgorithms)
//    useOldAlgorithm = true;

  crystal.omega=
    (crystal.bravais(1,0)*crystal.bravais(2,1)
     -crystal.bravais(2,0)*crystal.bravais(1,1))*crystal.bravais(0,2)+
    (crystal.bravais(2,0)*crystal.bravais(0,1)
     -crystal.bravais(0,0)*crystal.bravais(2,1))*crystal.bravais(1,2)+
    (crystal.bravais(0,0)*crystal.bravais(1,1)
     -crystal.bravais(1,0)*crystal.bravais(0,1))*crystal.bravais(2,2);
  crystal.omega=std::abs(crystal.omega);

// Testing the new VP construction:
  if(testingVorpol && (crystal.num_atoms < numSwitchAlgorithms) && !useOldAlgorithm)
  {
    std::vector<Real> atom_position_1(crystal.num_atoms);
    std::vector<Real> atom_position_2(crystal.num_atoms);
    std::vector<Real> atom_position_3(crystal.num_atoms);
    std::vector<Real> rad(crystal.num_atoms);
    for(int i=0; i<crystal.num_atoms; i++)
    {
      atom_position_1[i]=crystal.position(0,i);
      atom_position_2[i]=crystal.position(1,i);
      atom_position_3[i]=crystal.position(2,i);
      rad[i] = crystal.types[crystal.type[i]].rad;
    }

    for(int i=0; i<local.num_local; i++)
    {
      int my_atom = local.global_id[i]+1;
      int num_atoms=crystal.num_atoms;
      setup_boundary_(&my_atom, &num_atoms,
                      &atom_position_1[0],
                      &atom_position_2[0], &atom_position_3[0],
                      &crystal.bravais(0,0),
                      &crystal.bravais(0,1), &crystal.bravais(0,2),
                      &vplanesTest(0,0,i), &ipvp, &nvplanesTest[i], &rad[0]);
    }
  }

  if(useOldAlgorithm)
  {
    std::vector<Real> atom_position_1(crystal.num_atoms);
    std::vector<Real> atom_position_2(crystal.num_atoms);
    std::vector<Real> atom_position_3(crystal.num_atoms);
    std::vector<Real> rad(crystal.num_atoms);
    for(int i=0; i<crystal.num_atoms; i++)
    {
      atom_position_1[i]=crystal.position(0,i);
      atom_position_2[i]=crystal.position(1,i);
      atom_position_3[i]=crystal.position(2,i);
      rad[i] = crystal.types[crystal.type[i]].rad;
    }

    for(int i=0; i<local.num_local; i++)
    {
      int my_atom = local.global_id[i]+1;
      int num_atoms=crystal.num_atoms;
      setup_boundary_(&my_atom, &num_atoms,
                      &atom_position_1[0],
                      &atom_position_2[0], &atom_position_3[0],
                      &crystal.bravais(0,0),
                      &crystal.bravais(0,1), &crystal.bravais(0,2),
                      &vplanes(0,0,i), &ipvp, &nvplanes[i], &rad[0]);
    }

  } else {
    std::vector<Real> atom_position_1(maxClusterSize);
    std::vector<Real> atom_position_2(maxClusterSize);
    std::vector<Real> atom_position_3(maxClusterSize);
    std::vector<Real> rad(maxClusterSize);
    for(int i=0; i<local.num_local; i++)
    {
      int my_atom = 1;
      int num_atoms = local.atom[i].vpClusterGlobalIdx.size();
      for(int j=0; j<num_atoms; j++)
      {
        int idx = local.atom[i].vpClusterGlobalIdx[j];
        atom_position_1[j] = local.atom[i].vpClusterPos(0, j);
        atom_position_2[j] = local.atom[i].vpClusterPos(1, j);
        atom_position_3[j] = local.atom[i].vpClusterPos(2, j);
        rad[j] = crystal.types[crystal.type[idx]].rad;
      }

      setup_boundary_cluster_(&my_atom, &num_atoms,
                      &atom_position_1[0],
                      &atom_position_2[0], &atom_position_3[0],
                              &vplanes(0,0,i), &ipvp, &nvplanes[i], &rad[0]);
// Testing:
      if(testingVorpol && (crystal.num_atoms < numSwitchAlgorithms) && (lsms.global.iprint>=0))
      {
        printf("Testing Voronoi Boundery Planes [local = %d]:\n\n",i);
        printf("nvplanesTest = %d     nvplanes = %d\n", nvplanesTest[i], nvplanes[i]);
        for(int j=0; j<nvplanesTest[i]; j++)
          printf("vplanesTest(:,%d,%d) = %f %f %f\n",j,i,vplanesTest(0,j,i),vplanesTest(1,j,i),vplanesTest(2,j,i));
        for(int j=0; j<nvplanes[i]; j++)
          printf("vplanes(:,%d,%d) = %f %f %f\n",j,i,vplanes(0,j,i),vplanes(1,j,i),vplanes(2,j,i));
      }
    }
  }

  for(int i=0; i<local.num_local; i++)
  {
    int lmax=2*local.atom[i].lmax;
    local.atom[i].rInscribed=-1.0;
    local.atom[i].voronoi.rInscribedSphere=-1.0;
    local.atom[i].voronoi.wylm.resize((2*lmax+1)*(lmax+1),lsms.ngaussr,iprcrit-1);
    local.atom[i].voronoi.gwwylm.resize(lsms.ngaussr,iprcrit-1);
    local.atom[i].voronoi.grwylm.resize(lsms.ngaussr,iprcrit-1);
    int my_atom=local.global_id[i]+1;
    int num_atoms=crystal.num_atoms;
    /*
    setup_vorpol_(&my_atom,&num_atoms,
                  &atom_position_1[0],
                  &atom_position_2[0], &atom_position_3[0],
                  &crystal.bravais(0,0),
                  &lmax,&shc.clm[0],&lsms.ngaussq,&lsms.ngaussr,
                  &local.atom[i].voronoi.rInscribedSphere,&local.atom[i].voronoi.omegaInt,
                  local.atom[i].voronoi.dipint,&rad[0],
                  &ipvp,&ipnode,&ipcorn,&ipedge,&iprcrit,
                  &local.atom[i].voronoi.gwwylm(0,0),&local.atom[i].voronoi.grwylm(0,0),
                  &local.atom[i].voronoi.ncrit,&local.atom[i].voronoi.wylm(0,0,0),
                  &local.atom[i].rCircumscribed,
                  &lsms.global.iprint,lsms.global.istop,32);
    */

    setup_vorpol_vplane_(&vplanes(0,0,i), &nvplanes[i],
                         &lmax,&SphericalHarmonicsCoeficients::clm[0],&lsms.ngaussq,&lsms.ngaussr,
                         &local.atom[i].voronoi.rInscribedSphere,&local.atom[i].voronoi.omegaInt,
                         local.atom[i].voronoi.dipint,
                         &ipvp,&ipnode,&ipcorn,&ipedge,&iprcrit,
                         &local.atom[i].voronoi.gwwylm(0,0),&local.atom[i].voronoi.grwylm(0,0),
                         &local.atom[i].voronoi.ncrit,&local.atom[i].voronoi.wylm(0,0,0),
                         &local.atom[i].rCircumscribed,
                         &lsms.global.iprint,lsms.global.istop,32);

    local.atom[i].rInscribed=local.atom[i].voronoi.rInscribedSphere;
// set rmt according to value of fixRMT
    if(lsms.fixRMT==0)
    {
      local.atom[i].rmt=local.atom[i].rInscribed;
      local.atom[i].generateNewMesh = true;
    }
    // do we need to change rmt for mtasa==1?

    local.atom[i].omegaMT=sphereVolumeFactor*std::pow(local.atom[i].rmt,3);
    local.atom[i].omegaWS=local.atom[i].voronoi.omegaInt+local.atom[i].omegaMT;
    local.atom[i].rws=std::pow(local.atom[i].omegaWS/sphereVolumeFactor,1.0/3.0);

    switch(lsms.mtasa)
    {
    case 1:
      local.atom[i].rmt=local.atom[i].rws;
      local.atom[i].omegaMT=local.atom[i].omegaWS;
      local.atom[i].rInscribed=local.atom[i].rws;
      break;
    case 2:
      local.atom[i].rmt=local.atom[i].rws;
      local.atom[i].omegaMT=local.atom[i].omegaWS;
      local.atom[i].rInscribed=local.atom[i].voronoi.rInscribedSphere;
      break;
    default: // MT
      local.atom[i].rInscribed=local.atom[i].voronoi.rInscribedSphere;
    }
 }

  timeSetupVorpol = MPI_Wtime() - timeSetupVorpol;
  if(lsms.global.iprint >= 0)
    printf("time setupVorpol : %lf sec\n", timeSetupVorpol);
}


void calculateVolumes(LSMSCommunication &comm, LSMSSystemParameters &lsms, CrystalParameters &crystal, LocalTypeInfo &local)
{
  Real volumeTotal = 0.0;
  Real volumeMT = 0.0;
  const Real sphereVolumeFactor = 4.0 * M_PI / 3.0;

  for (int i=0; i<local.num_local; i++)
  {
    volumeTotal += local.atom[i].omegaWS * local.n_per_type[i];
    switch (lsms.mtasa)
    {
      case 1:
        volumeMT += local.atom[i].omegaMT;
        break;
      default:
        volumeMT += sphereVolumeFactor * std::pow(local.atom[i].rInscribed, 3);
    }
  }

  globalSum(comm, volumeTotal);
  globalSum(comm, volumeMT);

  lsms.volumeTotal = volumeTotal;
  lsms.volumeNorm = crystal.omega / volumeTotal;
  lsms.volumeInterstitial = volumeTotal - volumeMT;

  if(lsms.global.iprint >= 0)
  {
    std::printf("\n");
    std::printf("Total cell volume     = %20.13f\n", lsms.volumeTotal);
    std::printf("Volume renorm. factor = %20.13f\n", lsms.volumeNorm);
    std::printf("WS cell volume        = %20.13f\n", local.atom[0].omegaWS);
    std::printf("Interstitial volume   = %20.13f\n", lsms.volumeInterstitial);
    std::printf("WS sphere radius      = %20.13f\n", local.atom[0].rws);
    std::printf("Circumscribed radius  = %20.13f\n", local.atom[0].rCircumscribed);
  }
}

