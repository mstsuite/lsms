/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <chrono>
#include <cfenv>

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

#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/printf.h>
#include <hdf5.h>
#include <lua.hpp>


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
#include "num_digits.hpp"

#include "ChargeDensity/ChargeDensity.hpp"

#ifdef USE_NVTX
#include <nvToolsExt.h>
#endif


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

  auto lsmsStartTime = std::chrono::steady_clock::now();

  lua_State *L = luaL_newstate();
  luaL_openlibs(L);
  initLSMSLuaInterface(L);

  // feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

#ifdef USE_GPTL
  GPTLinitialize();
#endif
  initializeCommunication(comm);
  H5open();

  // set input file name (default 'i_lsms')
  strncpy(inputFileName, "i_lsms", 10);
  if (argc > 1)
    strncpy(inputFileName, argv[1], 120);

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
    printf("LSMS_3: Program started\n");
    printf("Using %d MPI processes\n", comm.size);
#ifdef _OPENMP
    printf("Using %d OpenMP threads\n", omp_get_max_threads());
#endif
    acceleratorPrint();
#ifdef LSMS_NO_COLLECTIVES
    printf("\nWARNING!!!\nCOLLECTIVE COMMUNICATION (ALLREDUCE etc.) ARE SKIPPED!\n");
    printf("THIS IS FOR TESTING ONLY!\nRESULTS WILL BE WRONG!!!\n\n");
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
    switch (lsms.mtasa)
    {
      case 1:
        printf("Performing Atomic Sphere Approximation (ASA) calculation\n");
        break;
      case 2:
        printf("Performing Atomic Sphere Approximation + Muffin-Tin (ASA-MT) calculation\n");
        break;
      default:
        printf("Performing Muffin-Tin (MT) calculation\n");
    }
    fflush(stdout);
  }

#ifdef LSMS_DEBUG
  MPI_Barrier(comm.comm);
#endif

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

#if defined(ACCELERATOR_CUDA_C) || defined(ACCELERATOR_HIP)
  deviceAtoms.resize(local.num_local);
#endif

  if(comm.rank == 0)
  {
    printf("set global ids.\n");
    fflush(stdout);
  }

#ifdef LSMS_DEBUG
  MPI_Barrier(comm.comm);
// write the distribution of atoms
  std::vector<int> atomsPerNode(comm.size);

#endif

  // set up exchange correlation functionals
  if (lsms.xcFunctional[0] == 1)  {
    lsms.exch_corr = std::make_shared<lsms::XCLibxc>(lsms.n_spin_pola, lsms.xcFunctional);
  }       // use libxc functional
  if (lsms.xcFunctional[0] == 2) {         // use new LSMS functional
    lsms.newFunctional.init(lsms.n_spin_pola, lsms.xcFunctional);
  }

  AngularMomentumIndices::init(2*crystal.maxlmax);
  SphericalHarmonicsCoeficients::init(2*crystal.maxlmax);

  GauntCoeficients::init(lsms);
  IFactors::init(lsms, crystal.maxlmax);

#if defined(ACCELERATOR_CUDA_C) || defined(ACCELERATOR_HIP)
  deviceConstants.allocate();
#endif

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

#ifdef LSMS_DEBUG
  MPI_Barrier(comm.comm);
#endif

// initialize the potential accelerators (GPU)
// we need to know the max. size of the kkr matrix to invert: lsms.n_spin_cant*local.maxNrmat()
// which is only available after building the LIZ

  acceleratorInitialize(lsms.n_spin_cant*local.maxNrmat(), lsms.global.GPUThreads);
  local.tmatStore.pinMemory();
#if defined(ACCELERATOR_CUBLAS) || defined(ACCELERATOR_LIBSCI) || defined(ACCELERATOR_CUDA_C) || defined(ACCELERATOR_HIP)
  // deviceStorage = allocateDStore();
  deviceStorage = new DeviceStorage;
#endif

  for (int i=0; i<local.num_local; i++)
    local.atom[i].pmat_m.resize(lsms.energyContour.groupSize());

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

//  initialAtomSetup(comm,lsms,crystal,local);

// the next line is a hack for initialization of potentials from scratch to work.


#ifdef LSMS_DEBUG
  if(lsms.global.iprint >= 0)
  {
    printf("Entering the Voronoi construction BEFORE loading the potentials.\n");
    fflush(stdout);
  }
  MPI_Barrier(comm.comm);
#endif

  /* if(lsms.pot_in_type < 0) */ setupVorpol(lsms, crystal, local);

#ifdef LSMS_DEBUG
  if(lsms.global.iprint >= 0)
  {
    printf("Entering the LOADING of the potentials.\n");
    fflush(stdout);
  }
  MPI_Barrier(comm.comm);
#endif

  double timeLoadPotential = MPI_Wtime();
  loadPotentials(comm, lsms, crystal, local);
  timeLoadPotential = MPI_Wtime() - timeLoadPotential;
  
  if (lsms.global.iprint >= 0)
  {
    printf("time loadPotential: %lf sec\n\n",timeLoadPotential);
    
    fprintf(stdout,"LIZ for atom 0 on this node\n");
    printLIZInfo(stdout, local.atom[0]);
    if(local.atom[0].forceZeroMoment) {
      fprintf(stdout,"\nMagnetic moment of atom 0 forced to be zero!\n\n");
    }
  }

  // Read evec file if in is define
  if ( lsms.infoEvecFileIn[0]!=0) {
    readInfoEvec(comm,
                 lsms,
                 crystal,
                 local,
                 lsms.infoEvecFileIn);
    if (lsms.global.iprint >= 0)
    {
      fprintf(stdout,"Evec are read from: %s\n", lsms.infoEvecFileIn);
    }
  }


  if ( !alloyDesc.empty() )
  {
    if(lsms.global.iprint >= 0)
    {
      printf("Entering the LOADING of the alloy banks.\n");
      fflush(stdout);
    }
    loadAlloyBank(comm,lsms,alloyDesc,alloyBank);
  }

// for testing purposes:
//  std::vector<Matrix<Real> > vrs;
//  vrs.resize(local.num_local);
//  for(int i=0; i<local.num_local; i++) vrs[i]=local.atom[i].vr;
// -------------------------------------

#ifdef LSMS_DEBUG
  if(lsms.global.iprint >= 0)
  {
    printf("Entering the Voronoi construction AFTER loading the potentials.\n");
    fflush(stdout);
  }
  MPI_Barrier(comm.comm);
#endif

  setupVorpol(lsms, crystal, local);

#ifdef LSMS_DEBUG
  MPI_Barrier(comm.comm);
#endif

// Generate new grids after new rmt is defined
  for (int i=0; i<local.num_local; i++)
  {
    if(local.atom[i].generateNewMesh)
      interpolatePotential(lsms, local.atom[i]);
  }

  double timeCalculateVolumes = MPI_Wtime();
  calculateVolumes(comm, lsms, crystal, local);
  timeCalculateVolumes = MPI_Wtime() - timeCalculateVolumes;
  if (lsms.global.iprint >= 0)
  {
    printf("time calculateVolumes: %lf sec\n\n",timeCalculateVolumes);
  }
  
//  loadPotentials(comm,lsms,crystal,local);
  if (lsms.pot_in_type == -1) {

    std::vector<Real> qsub(crystal.num_types, 0.0);

    Array3d<Real> rhoTemp;
    rhoTemp.resize(lsms.global.iprpts + 1, 2, local.num_local);
    rhoTemp = 0.0;
    
    calculatePotential(comm, lsms, local, crystal, qsub, rhoTemp, 0);

    // Initialize potentials and charge densities
    lsms::copyChargesAndPotential(lsms, local);

    // Check charge density after mixing
    lsms::checkRadialChargeDensity(lsms, local);

    if (comm.rank == 0) {
      fmt::printf("Initial MTZ: %20.9f\n", lsms.vmt);
    }

  }

// initialize Mixing
  double timeSetupMixing = MPI_Wtime();
  Mixing *mixing;
  setupMixing(mix, mixing, lsms.global.iprint);
  timeSetupMixing = MPI_Wtime() - timeSetupMixing;
  if (lsms.global.iprint >= 0)
  {
    printf("time setupMixing: %lf sec\n\n",timeSetupMixing);
  }

  double timeCalculateMadelungMatrix = MPI_Wtime();
// need to calculate madelung matrices
//#define LEGACY_MONOPOLE
#ifdef LEGACY_MONOPOLE
  calculateMadelungMatrices(lsms, crystal, local);
#else
  lsms::calculateMultiMadelungMatrices(lsms, crystal, local);
#endif
  timeCalculateMadelungMatrix = MPI_Wtime() - timeCalculateMadelungMatrix;
  if (lsms.global.iprint >= 0)
  {
    printf("time calculateMultiMadelungMatrices: %lf sec\n\n",timeCalculateMadelungMatrix);
  }

  if (lsms.global.iprint >= 1)
  {
    printLocalTypeInfo(stdout, local);
  }

  lsms::calculateCoreStates(comm, lsms, local);
  if (lsms.global.iprint >= 0)
  {
    printf("Finished calculateCoreStates(...)\n");
    fflush(stdout);
  }

// check that vrs have not changed ...
//  bool vr_check=false;
//  for(int i=0; i<local.num_local; i++)
//  {
//    vr_check=true;
//    for(int j=0; j<vrs[i].n_row();j++)
//      for(int k=0; k<vrs[i].n_col(); k++)
//        vr_check=vr_check && (vrs[i](j,k)==local.atom[i].vr(j,k));
//    if(!vr_check)
//      printf("Potential %d has changed!!\n",i);
//  }
//  printf("Potentials have been checked\n");
// --------------------------------------------

// meis: Test if Potential for atoms 0 and 1 are the same
/*
  if(local.num_local>1)
  {
    bool vr_check=true;
    for(int j=0; j<local.atom[0].vr.n_row();j++)
      for(int k=0; k<local.atom[0].vr.n_col(); k++)
        vr_check=vr_check && (local.atom[0].vr(j,k)==local.atom[1].vr(j,k));
    if(!vr_check)
      printf("Potentials 0 and 1 are different!!\n");
    printf("Potentials have been checked\n");
  }
*/

  if (lsms.n_spin_cant > 1)
  {
    for (int i=0; i<local.num_local; i++)
      local.atom[i].get_b_basis();
  }
  else
  {
    for (int i=0; i<local.num_local; i++)
      local.atom[i].reset_b_basis();
  }

  double timeMixingPrepare = MPI_Wtime();
  mixing -> prepare(comm, lsms, local.atom);
  timeMixingPrepare = MPI_Wtime() - timeMixingPrepare;
  if (lsms.global.iprint >= 0)
  {
    printf("time Mixing->Prepare: %lf sec\n\n",timeMixingPrepare);
  }

#ifdef USE_PAPI
  #define NUM_PAPI_EVENTS 2
  int hw_counters = PAPI_num_counters();
  if (hw_counters > NUM_PAPI_EVENTS)
    hw_counters = NUM_PAPI_EVENTS;
  int papi_events[NUM_PAPI_EVENTS]; 
  char *papi_event_name[] = {"PAPI_TOT_INS", "PAPI_FP_OPS"};
  // get events from names:
  for (int i=0; i<NUM_PAPI_EVENTS; i++)
  {
    if (PAPI_event_name_to_code(papi_event_name[i], &papi_events[i]) != PAPI_OK)
      if (hw_counters > i)
        hw_counters = i;
  }
  long long papi_values[NUM_PAPI_EVENTS+4];
  if (hw_counters > NUM_PAPI_EVENTS)
    hw_counters = NUM_PAPI_EVENTS;
  long long papi_real_cyc_0 = PAPI_get_real_cyc();
  long long papi_real_usec_0 = PAPI_get_real_usec();
  long long papi_virt_cyc_0 = PAPI_get_virt_cyc();
  long long papi_virt_usec_0 = PAPI_get_virt_usec();
  PAPI_start_counters(papi_events, hw_counters);
#endif

  auto lsmsEndInitTime = std::chrono::steady_clock::now();
  
// -----------------------------------------------------------------------------
//                                 MAIN SCF LOOP
// -----------------------------------------------------------------------------

  bool converged = false;
  bool energyConverged = true;
  Real oldTotalEnergy = lsms.totalEnergy;

  if (lsms.global.iprint >= 0)
  {
    printf("Total number of iterations:%d\n", lsms.nscf);
    fflush(stdout);
  }
    
  double timeScfLoop = MPI_Wtime();
  double timeCalcChemPot = 0.0;
  double timeCalcPotentialsAndMixing = 0.0;

  int iterationStart = 0;
  int potentialWriteCounter = 0;

  FILE *kFile = nullptr;
  if (comm.rank == 0)
  {
    iterationStart = readNextIterationNumber("k.out");
    kFile = fopen("k.out","a");
  }

  int iteration;
#ifdef USE_NVTX
  nvtxRangePushA("SCFLoop");
#endif
  for (iteration=0; iteration<lsms.nscf && !(converged && energyConverged); iteration++)
  {
    if (lsms.global.iprint >= -1 && comm.rank == 0)
      printf("SCF iteration %d:\n", iteration);

    oldTotalEnergy = lsms.totalEnergy;

    // Calculate band energy
    energyContourIntegration(comm, lsms, local);
    double dTimeCCP = MPI_Wtime();
    // if(!lsms.global.checkIstop("buildKKRMatrix"))

    // Calculate chemical potential 
    lsms::calculateChemPot(comm, lsms, local, eband);
    dTimeCCP = MPI_Wtime() - dTimeCCP;
    timeCalcChemPot += dTimeCCP;

    // Calculate magnetic moments for each site and check if spin has flipped
    calculateEvec(lsms, local);
    // mixEvec(lsms, local, 0.0);
    mixing -> updateMoments(comm, lsms, local.atom);
    for (int i=0; i<local.num_local; i++)
    {
      if(!mix.quantity[MixingParameters::moment_direction])
        local.atom[i].newConstraint();

      local.atom[i].evec[0] = local.atom[i].evecNew[0];
      local.atom[i].evec[1] = local.atom[i].evecNew[1];
      local.atom[i].evec[2] = local.atom[i].evecNew[2];

      checkIfSpinHasFlipped(lsms, local.atom[i]);
    }

    double dTimePM = MPI_Wtime();
    // Calculate charge densities, potentials, and total energy
    calculateAllLocalChargeDensities(lsms, local);
    calculateChargesPotential(comm, lsms, local, crystal, 0);

    // FILE *pf = fopen("vr_test_1.dat","w");
    // printAtomPotential(pf, local.atom[0]);
    // fclose(pf);

    checkAllLocalCharges(lsms, local);
    calculateTotalEnergy(comm, lsms, local, crystal);

    // pf = fopen("vr_test_2.dat","w");
    // printAtomPotential(pf, local.atom[0]);
    // fclose(pf);

    // Calculate charge density rms
    calculateLocalQrms(lsms, local);

    // Mix charge density
    mixing -> updateChargeDensity(comm, lsms, local.atom);
    dTimePM = MPI_Wtime() - dTimePM;
    timeCalcPotentialsAndMixing += dTimePM;

    // pf = fopen("vr_test_3.dat","w");
    // printAtomPotential(pf, local.atom[0]);
    // fclose(pf);

    // Recalculate core states
    // - swap core state energies for different spin channels first if spin has flipped
    //   (from LSMS 1: lsms_main.f:2101-2116)
    for (int i=0; i<local.num_local; i++) {
      if (local.atom[i].spinFlipped)
      {
        checkIfSpinHasFlipped(lsms, local.atom[i]);
        if (!local.atom[i].spinFlipped)
          swapCoreStateEnergies(local.atom[i]);
      }
    }
    lsms::calculateCoreStates(comm, lsms, local);

    // pf = fopen("vr_test_4.dat","w");
    // printAtomPotential(pf, local.atom[0]);
    // fclose(pf);

    dTimePM = MPI_Wtime();
    // If charge is mixed, recalculate potential and mix (need a flag for this from input)
    calculateChargesPotential(comm, lsms, local, crystal, 1);

    // pf = fopen("vr_test_5.dat","w");
    // printAtomPotential(pf, local.atom[0]);
    // fclose(pf);

    mixing -> updatePotential(comm, lsms, local.atom);

    // pf = fopen("vr_test_6.dat","w");
    // printAtomPotential(pf, local.atom[0]);
    // fclose(pf);

    // exitLSMS(comm, 1);

    dTimePM = MPI_Wtime() - dTimePM;
    timeCalcPotentialsAndMixing += dTimePM;

    Real rms = 0.0;
    if(lsms.n_spin_pola > 1)
    {
      // rms = 0.5 * (local.qrms[0] + local.qrms[1]);
      rms = 0.0;
      for(int i=0; i<local.num_local; i++)
        rms = std::max(rms, 0.5*(local.atom[i].qrms[0]+local.atom[i].qrms[1]));
      globalMax(comm, rms);
    } else {
      // rms = local.qrms[0];
      rms = 0.0;
      for(int i=0; i<local.num_local; i++)
        rms = std::max(rms, local.atom[i].qrms[0]);
      globalMax(comm, rms);
    }

// check for convergence
    converged = rms < lsms.rmsTolerance;
    /*
    converged = true;
    for (int i=0; i<local.num_local; i++)
    {
      converged = converged
                && (0.5*(local.atom[i].qrms[0]+local.atom[i].qrms[1])<lsms.rmsTolerance);
    }
    globalAnd(comm, converged);
    */
    if(lsms.energyTolerance > 0)
      energyConverged = std::abs((lsms.totalEnergy - oldTotalEnergy)/lsms.totalEnergy) < lsms.energyTolerance;
    else
      energyConverged = true;

    if (comm.rank == 0)
    {

      int gap_size = 12;
      int size = -1;
      size = std::max(size, num_digits(static_cast<int> (eband)));
      size = std::max(size, num_digits(static_cast<int> (lsms.totalEnergy)));
      size = std::max(size, num_digits(static_cast<int> (lsms.chempot)));
      size = std::max(size, num_digits(static_cast<int> (lsms.vmt)));

      gap_size -= size;
      gap_size = std::max(gap_size, 2);
      size += 10;

      std::printf("MTZ          = %*.9f Ry\n", size, lsms.vmt);
      std::printf("Band Energy  = %*.9f Ry %*s Fermi Energy = %12.9f Ry\n", size, eband,
                  gap_size, "", lsms.chempot);
      std::printf("Total Energy = %*.9f Ry\n", size, lsms.totalEnergy);
      std::printf("RMS = %lg\n",rms);
      if(lsms.global.iprint > 0)
      {
        printf("  qrms[0] = %lg   qrms[1] = %lg\n",local.qrms[0], local.qrms[1]);
        printf("  local.atom[i]:\n");
        for (int i=0; i<local.num_local; i++)
        {
          printf("  %d : qrms[0] = %lg   qrms[1] = %lg\n",i,local.atom[i].qrms[0], local.atom[i].qrms[1]);
          printf("  %d : vrms[0] = %lg   vrms[1] = %lg\n",i,local.atom[i].vrms[0], local.atom[i].vrms[1]);
        }
      }
    }

    if (kFile != nullptr)
    {
      fprintf(kFile,"%4d %20.12lf %12.6lf %12.6lf  %14.10lf\n",
              iterationStart+iteration, lsms.totalEnergy, lsms.chempot, local.atom[0].mtotws, rms);
      fflush(kFile);
    }

    // Recalculate core states for new potential if we are performing scf calculations
    lsms::calculateCoreStates(comm, lsms, local);

    // Periodically write the new potential for scf calculations 
    potentialWriteCounter++;
    if ((lsms.pot_out_type >= 0 && potentialWriteCounter >= lsms.writeSteps)
        || converged)
    {
      if (comm.rank == 0) std::cout << "Writing new potentials and restart file.\n";
      writePotentials(comm, lsms, crystal, local);
      potentialWriteCounter = 0;
      if (comm.rank == 0)
      {
        writeRestart("i_lsms.restart", lsms, crystal, mix, potentialShifter, alloyDesc);
      }
    }

  }
#ifdef USE_NVTX
  nvtxRangePop();
#endif
  timeScfLoop = MPI_Wtime() - timeScfLoop;

  writeInfoEvec(comm, lsms, crystal, local, eband, lsms.infoEvecFileOut);
  if(lsms.localAtomDataFile[0]!=0)
    writeLocalAtomData(comm, lsms, crystal, local, eband, lsms.localAtomDataFile);

  if (kFile != nullptr)
    fclose(kFile);

  /**
   * Total energy calculation of all contributions
   */

  lsms::DFTEnergy dft_energy;
  calculateTotalEnergy(comm, lsms, local, crystal, dft_energy);

  if (comm.rank == 0)
  {
    lsms::print_dft_energy( dft_energy);
  }

// -----------------------------------------------------------------------------

#ifdef USE_PAPI
  PAPI_stop_counters(papi_values,hw_counters);
  papi_values[hw_counters  ] = PAPI_get_real_cyc()-papi_real_cyc_0;
  papi_values[hw_counters+1] = PAPI_get_real_usec()-papi_real_usec_0;
  papi_values[hw_counters+2] = PAPI_get_virt_cyc()-papi_virt_cyc_0;
  papi_values[hw_counters+3] = PAPI_get_virt_usec()-papi_virt_usec_0;
  long long accumulated_counters[NUM_PAPI_EVENTS+4];
  MPI_Reduce(papi_values,accumulated_counters,hw_counters+4,
             MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
  if (comm.rank == 0)
  {
    for (int i=0; i<hw_counters; i++)
      std::cout<<"Accumulated: "<<(papi_event_name[i])<<" = "<<(accumulated_counters[i])<<"\n";
    std::cout<<"PAPI accumulated real cycles : "<<(accumulated_counters[hw_counters])<<"\n";
    std::cout<<"PAPI accumulated user cycles : "<<(accumulated_counters[hw_counters+2])<<"\n";
    double gflops_papi = ((double)accumulated_counters[1])/
      (1000.0*(double)papi_values[hw_counters+1]);
    double gflops_hw_double = ((double)accumulated_counters[2])/
      (1000.0*(double)papi_values[hw_counters+1]);
    double gflops_hw_single = ((double)accumulated_counters[3])/
      (1000.0*(double)papi_values[hw_counters+1]);
    double gips = ((double)accumulated_counters[0])/(1000.0*(double)papi_values[hw_counters+1]);
    std::cout<<"PAPI_FP_OPS real GFLOP/s : "<<(gflops_papi)<<"\n";
    std::cout<<"PAPI hw double real GFLOP/s : "<<(gflops_hw_double)<<"\n";
    std::cout<<"PAPI hw single real GFLOP/s : "<<(gflops_hw_single)<<"\n";
    std::cout<<"PAPI real GINST/s : "<<(gips)<<"\n";
    std::cout<<"Time (s) : " << (double)papi_values[hw_counters+1] << "\n";
  }
#endif

  if (lsms.pot_out_type >= 0)
  {
    if (comm.rank == 0) std::cout << "Writing new potentials.\n";
    writePotentials(comm, lsms, crystal, local);
    if (comm.rank == 0)
    {
      std::cout << "Writing restart file.\n";
      writeRestart("i_lsms.restart", lsms, crystal, mix, potentialShifter, alloyDesc);
    }
  }

  double fomScale = calculateFomScaleDouble(comm, local);

  auto lsmsEndTime = std::chrono::steady_clock::now();
  std::chrono::duration<double> lsmsRuntime = lsmsEndTime - lsmsStartTime;
  std::chrono::duration<double> lsmsInitTime = lsmsEndInitTime - lsmsStartTime;
  
  if (comm.rank == 0)
  {
    int size = -1;
    size = std::max(size, num_digits(static_cast<int> (eband)));
    size = std::max(size, num_digits(static_cast<int> (lsms.chempot)));
    size = std::max(size, num_digits(static_cast<int> (lsms.totalEnergy)));
    size += 17;

    printf("Band Energy  = %*.15f Ry\n", size, eband);
    printf("Fermi Energy = %*.15f Ry\n", size, lsms.chempot);
    printf("Total Energy = %*.15f Ry\n", size, lsms.totalEnergy);
    printf("\nTimings:\n========\n");
    printf("LSMS Runtime = %lf sec\n", lsmsRuntime.count());
    printf("LSMS Initialization Time = %lf sec\n", lsmsInitTime.count());
    printf("timeScfLoop[rank==0] = %lf sec\n", timeScfLoop);
    printf("     number of iteration:%d\n",iteration);
    printf("timeScfLoop/iteration = %lf sec\n", timeScfLoop / (double)iteration);
    // printf(".../lsms.nscf = %lf sec\n", timeScfLoop / (double)lsms.nscf);
    printf("timeCalcChemPot[rank==0]/iteration = %lf sec\n", timeCalcChemPot / (double)iteration);
    // printf("timeCalcChemPot[rank==0]/lsms.nscf = %lf sec\n", timeCalcChemPot / (double)lsms.nscf);
    printf("timeCalcPotentialsAndMixing[rank==0]/iteration = %lf sec\n",
            timeCalcPotentialsAndMixing / (double)iteration);
    printf("timeBuildLIZandCommList[rank==0]: %lf sec\n",
           timeBuildLIZandCommList);
    // fom = [ \sum_#atoms (LIZ * (lmax+1)^2)^3 ] / time per iteration
    //     = [ \sum_#atoms (LIZ * (lmax+1)^2)^3 ] * lsms.nscf / timeScfLoop
    // fom_e = fom * energy contour points
    // fomScale = \sum_#atoms (LIZ * (lmax+1)^2)^3
    // energyContourPoints
    long energyContourPoints = 1;
    if(lsms.energyContour.grid==2)
    {
      energyContourPoints = lsms.energyContour.npts+1;
    }
    printf("FOM Scale = %lf\n",(double)fomScale);
    printf("Energy Contour Points = %ld\n",energyContourPoints);
    printf("FOM / energyContourPoint = %lg/sec\n",fomScale * (double)iteration / timeScfLoop);
    // printf("FOM = %lg/sec\n",fomScale * (double)lsms.nscf / timeScfLoop);
    printf("FOM = %lg/sec\n",
            (double)energyContourPoints * (double)fomScale * (double)iteration / timeScfLoop);
    //         (double)energyContourPoints * (double)fomScale * (double)lsms.nscf / timeScfLoop);
  }


  local.tmatStore.unpinMemory();


#if defined(ACCELERATOR_CUBLAS) || defined(ACCELERATOR_LIBSCI) || defined(ACCELERATOR_CUDA_C) || defined(ACCELERATOR_HIP)
  // freeDStore(deviceStorage);
  delete deviceStorage;
#endif

  acceleratorFinalize();
  delete mixing;

#ifdef USE_GPTL
  GPTLpr(comm.rank);
#endif

  H5close();
  finalizeCommunication();
  lua_close(L);
  return EXIT_SUCCESS;
}
