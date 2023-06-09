/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

#include <fenv.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// #define USE_PAPI 1
#ifdef USE_PAPI
#include <papi.h>
#endif

#ifdef USE_GPTL
#include "gptl.h"
#endif

#include <hdf5.h>

#include "lua.hpp"


#include "LuaInterface/LuaInterface.hpp"
#include "SystemParameters.hpp"
#include "PotentialIO.hpp"
#include "Communication/distributeAtoms.hpp"
#include "Communication/LSMSCommunication.hpp"
#include "Core/calculateCoreStates.hpp"
#include "Misc/Indices.hpp"
#include "Misc/Coeficients.hpp"
#include "Madelung/Madelung.hpp"
#include "MultipoleMadelung/calculateMultipoleMadelung.hpp"
#include "VORPOL/VORPOL.hpp"
#include "energyContourIntegration.hpp"
#include "Accelerator/Accelerator.hpp"
#include "calculateChemPot.hpp"
#include "calculateDensities.hpp"
#include "mixing.hpp"
#include "calculateEvec.hpp"
#include "Potential/calculateChargesPotential.hpp"
#include "Potential/interpolatePotential.hpp"
#include "Potential/PotentialShifter.hpp"
#include "TotalEnergy/calculateTotalEnergy.hpp"
#include "SingleSite/checkAntiFerromagneticStatus.hpp"
#include "VORPOL/setupVorpol.hpp"
#include "Misc/readLastLine.hpp"
#include "Potential/XCBase.hpp"
#include "Potential/XCLibxc.hpp"

#include "buildLIZandCommLists.hpp"
#include "writeInfoEvec.hpp"
#include "write_restart.hpp"
#include "mixing_params.hpp"
#include "read_input.hpp"

#ifdef USE_NVTX
#include <nvToolsExt.h>
#endif

SphericalHarmonicsCoeficients sphericalHarmonicsCoeficients;
GauntCoeficients gauntCoeficients;
IFactors iFactors;

#if defined(ACCELERATOR_CUBLAS) || defined(ACCELERATOR_LIBSCI) || defined(ACCELERATOR_CUDA_C) || defined(ACCELERATOR_HIP)
#include "Accelerator/DeviceStorage.hpp"
// void * deviceStorage;
DeviceStorage *deviceStorage;
DeviceConstants deviceConstants;
#endif
// std::vector<void *> deviceConstants;
// std::vector<void *> deviceStorage;


/*
 * Need portablew way to enable FP exceptions!
 *
static int
feenableexcept (unsigned int excepts)
{
  static fenv_t fenv;
  unsigned int new_excepts = excepts & FE_ALL_EXCEPT,
               old_excepts;  // previous masks

  if ( fegetenv (&fenv) ) return -1;
  old_excepts = fenv.__control & FE_ALL_EXCEPT;

  // unmask
  fenv.__control &= ~new_excepts;
  fenv.__mxcsr   &= ~(new_excepts << 7);

  return ( fesetenv (&fenv) ? -1 : old_excepts );
}
*/ 


int main(int argc, char *argv[])
{
  LSMSSystemParameters lsms;
  LSMSCommunication comm;
  CrystalParameters crystal;
  LocalTypeInfo local;
  MixingParameters mix;
  PotentialShifter potentialShifter;
  AlloyMixingDesc alloyDesc;
  AlloyAtomBank alloyBank;

  char inputFileName[128];

  Real eband;

  lua_State *L = luaL_newstate();
  luaL_openlibs(L);
  initLSMSLuaInterface(L);

  initializeCommunication(comm);

  // simulate number of mpi ranks:
  int forceNumberOfMPIRanks = comm.size;

  // set input file name (default 'i_lsms')
  strncpy(inputFileName, "i_lsms", 10);
  if (argc > 1)
  {
    forceNumberOfMPIRanks = atoi(argv[1]);
    if (argc > 2)
      strncpy(inputFileName, argv[2], 120);
  }

  comm.size = forceNumberOfMPIRanks;

  lsms.global.iprpts = 1051;
  lsms.global.ipcore = 30;
  lsms.global.setIstop("main");
  lsms.global.iprint = 0;
  lsms.global.default_iprint = -1;
  lsms.global.print_node = 0;
  lsms.ngaussr = 10;
  lsms.ngaussq = 40;
  lsms.vSpinShiftFlag = 0;
#ifdef _OPENMP
  lsms.global.GPUThreads = std::min(12, omp_get_max_threads());
#else
  lsms.global.GPUThreads = 1;
#endif
  if (comm.rank == 0)
    lsms.global.iprint = 0;

  if (comm.rank == 0)
  {
    printf("LIZ and Voronoi Tewst: Program started\n");
    printf("Using %d MPI processes\n", comm.size);
#ifdef _OPENMP
    printf("Using %d OpenMP threads\n", omp_get_max_threads());
#endif

    printf("Reading input file '%s'\n", inputFileName);
    fflush(stdout);

    if (luaL_loadfile(L, inputFileName) || lua_pcall(L,0,0,0))
    {
      fprintf(stderr, "!! Cannot run input file!!\n");
      exit(1);
    }

    printf("Loaded input file!\n");
    fflush(stdout);

    if (readInput(L, lsms, crystal, mix, potentialShifter, alloyDesc))
    {
      fprintf(stderr, "!! Something wrong in input file!!\n");
      exit(1);
    }

    printf("System information:\n");
    printf("===================\n");
    printf("Number of atoms        : %10d\n", crystal.num_atoms);
    printf("Number of atomic types : %10d\n", crystal.num_types);

    fflush(stdout);
  }


  communicateParameters(comm, lsms, crystal, mix, alloyDesc);
  if (comm.rank != lsms.global.print_node)
    lsms.global.iprint = lsms.global.default_iprint;
  // printf("maxlmax=%d\n",lsms.maxlmax);
  if(comm.rank == 0)
  {
    printf("communicated Parameters.\n");
    fflush(stdout);
  }

  local.setNumLocal(distributeTypes(crystal, comm));
  local.setGlobalId(comm.rank, crystal);
     
  if(comm.rank == 0)
  {
    printf("set global ids.\n");
    fflush(stdout);
  }

  AngularMomentumIndices::init(2*crystal.maxlmax);
  SphericalHarmonicsCoeficients::init(2*crystal.maxlmax);

  GauntCoeficients::init(lsms);
  IFactors::init(lsms, crystal.maxlmax);

  double timeBuildLIZandCommList = MPI_Wtime();
  if (lsms.global.iprint >= 0)
  {
    printf("building the LIZ and Communication lists [buildLIZandCommLists]\n");
    fflush(stdout);
  }
  buildLIZandCommLists(comm, lsms, crystal, local);
  timeBuildLIZandCommList = MPI_Wtime() - timeBuildLIZandCommList;
  if (lsms.global.iprint >= 0)
  {
    printf("time for buildLIZandCommLists [num_local=%d]: %lf sec\n",
           local.num_local, timeBuildLIZandCommList);
    fflush(stdout);
  }

// set maximal number of radial grid points and core states if reading from bigcell file
  local.setMaxPts(lsms.global.iprpts);
  local.setMaxCore(lsms.global.ipcore);

  if (lsms.global.iprint >= 0) printLSMSGlobals(stdout, lsms);
  if (lsms.global.iprint >= 0) printLSMSSystemParameters(stdout, lsms);
  if (lsms.global.iprint >= 1) printCrystalParameters(stdout, crystal);
  if (lsms.global.iprint >= 0) printAlloyParameters(stdout, alloyDesc);

  if (lsms.global.iprint >= 1)
  {
    printCommunicationInfo(stdout, comm);
  }
  fflush(stdout);


  double timeVORPOL = MPI_Wtime();
  setupVorpol(lsms, crystal, local);
  timeVORPOL = MPI_Wtime() - timeVORPOL;
  if (lsms.global.iprint >= 0)
  {
    printf("time for setupVorpol [num_local=%d]: %lf sec\n",
           local.num_local, timeVORPOL);
    fflush(stdout);
  }

  double timeCalculateMadelungMatrix = MPI_Wtime();
// need to calculate madelung matrices
// #ifdef LEGACY_MONOPOLE
  if (lsms.global.iprint >= 0)
  {
    printf("Legacy Monopole Madelung Matrix\n");
  }
  calculateMadelungMatrices(lsms, crystal, local);
/*
#else
  if (lsms.global.iprint >= 0)
  {
    printf("New Multipole Madelung Matrix Routine\n");
  }
  calculateMultiMadelungMatrices(lsms, crystal, local);
#endif
*/
  timeCalculateMadelungMatrix = MPI_Wtime() - timeCalculateMadelungMatrix;
  if (lsms.global.iprint >= 0)
  {
    printf("time calculateMadelungMatrices: %lf sec\n\n",timeCalculateMadelungMatrix);
  }

  timeCalculateMadelungMatrix = MPI_Wtime();
// need to calculate madelung matrices
/*
#ifdef LEGACY_MONOPOLE
  if (lsms.global.iprint >= 0)
  {
    printf("Legacy Monopole Madelung Matrix\n");
  }
  calculateMadelungMatrices(lsms, crystal, local);
#else
*/
  if (lsms.global.iprint >= 0)
  {
    printf("New Multipole Madelung Matrix Routine\n");
  }
  lsms::calculateMultiMadelungMatrices(lsms, crystal, local);
// #endif
  timeCalculateMadelungMatrix = MPI_Wtime() - timeCalculateMadelungMatrix;
  if (lsms.global.iprint >= 0)
  {
    printf("time calculateMadelungMatrices: %lf sec\n\n",timeCalculateMadelungMatrix);
  }
  
  finalizeCommunication();
  lua_close(L);
  return 0;
}
