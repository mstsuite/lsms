
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

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

#include "Accelerator/Accelerator.hpp"
#include "ChargeDensity.hpp"
#include "ChargePlotter.hpp"
#include "Communication/LSMSCommunication.hpp"
#include "Communication/distributeAtoms.hpp"
#include "Core/calculateCoreStates.hpp"
#include "LuaInterface/LuaInterface.hpp"
#include "Main/setupCell.hpp"
#include "Misc/Coeficients.hpp"
#include "Misc/Indices.hpp"
#include "Misc/readLastLine.hpp"
#include "Mixer.hpp"
#include "MixingParameter.hpp"
#include "MultipoleMadelung/calculateMultipoleMadelung.hpp"
#include "Potential/PotentialShifter.hpp"
#include "Potential/XCBase.hpp"
#include "Potential/XCLibxc.hpp"
#include "Potential/interpolatePotential.hpp"
#include "PotentialIO.hpp"
#include "SystemParameters.hpp"
#include "TotalEnergy/calculateTotalEnergy.hpp"
#include "VORPOL/setupVorpol.hpp"
#include "buildLIZandCommLists.hpp"
#include "calculateChemPot.hpp"
#include "energyContourIntegration.hpp"
#include "newPotential.hpp"
#include "num_digits.hpp"
#include "read_input.hpp"
#include "updateChargePotential.hpp"
#include "writeInfoEvec.hpp"
#include "write_restart.hpp"
#include "Madelung.hpp"

#ifdef USE_NVTX
#include <nvToolsExt.h>
#endif

#if defined(ACCELERATOR_CUBLAS) || defined(ACCELERATOR_LIBSCI) || \
    defined(ACCELERATOR_CUDA_C) || defined(ACCELERATOR_HIP)
#include "Accelerator/DeviceStorage.hpp"
DeviceStorage *deviceStorage;
DeviceConstants deviceConstants;
#endif

int main(int argc, char *argv[]) {
  LSMSSystemParameters lsms;
  LSMSCommunication comm;
  CrystalParameters crystal;
  LocalTypeInfo local;
  lsms::MixingParameterPack mix;
  PotentialShifter potentialShifter;
  AlloyMixingDesc alloyDesc;
  AlloyAtomBank alloyBank;

  char inputFileName[128];

  Real eband;

  auto lsmsStartTime = std::chrono::steady_clock::now();

  lua_State *L = luaL_newstate();
  luaL_openlibs(L);

#ifdef USE_GPTL
  GPTLinitialize();
#endif
  initializeCommunication(comm);
  H5open();

  // set input file name (default 'i_lsms')
  strncpy(inputFileName, "i_lsms", 10);
  if (argc > 1) strncpy(inputFileName, argv[1], 120);

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
  if (comm.rank == 0) lsms.global.iprint = 0;

  if (comm.rank == 0) {
    fmt::printf("LSMS_3: Program started\n");
    fmt::printf("Using %d MPI processes\n", comm.size);
#ifdef _OPENMP
    fmt::printf("Using %d OpenMP threads\n", omp_get_max_threads());
#endif
    acceleratorPrint();
#ifdef LSMS_NO_COLLECTIVES
    fmt::printf(
        "\nWARNING!!!\nCOLLECTIVE COMMUNICATION (ALLREDUCE etc.) ARE "
        "SKIPPED!\n");
    fmt::printf("THIS IS FOR TESTING ONLY!\nRESULTS WILL BE WRONG!!!\n\n");
#endif
    fmt::printf("Reading input file '%s'\n", inputFileName);
    fflush(stdout);

    if (luaL_loadfile(L, inputFileName) || lua_pcall(L, 0, 0, 0)) {
      fprintf(stderr, "!! Cannot run input file!!\n");
      exit(1);
    }
    fmt::printf("Loaded input file!\n");
    fflush(stdout);

    if (readInput(L, lsms, crystal, mix, potentialShifter, alloyDesc)) {
      fprintf(stderr, "!! Something wrong in input file!!\n");
      exit(1);
    }

    // Fix to ASA
    lsms.mtasa = 1;

    fmt::printf("System information:\n");
    fmt::printf("Number of atoms        : %10d\n", crystal.num_atoms);
    fmt::printf("Number of atomic types : %10d\n", crystal.num_types);
    fmt::printf("Performing Atomic Sphere Approximation (ASA) calculation\n");
    fflush(stdout);
  }

#ifdef LSMS_DEBUG
  MPI_Barrier(comm.comm);
#endif

  communicateParameters(comm, lsms, crystal, mix, alloyDesc);
  if (comm.rank != lsms.global.print_node)
    lsms.global.iprint = lsms.global.default_iprint;

  if (comm.rank == 0) {
    fmt::printf("communicated Parameters.\n");
    fflush(stdout);
  }

  // Set the number of local types
  int nr_local_types = distributeTypes(crystal, comm);
  local.setNumLocal(nr_local_types);
  local.setGlobalId(comm.rank, crystal);

#if defined(ACCELERATOR_CUDA_C) || defined(ACCELERATOR_HIP)
  deviceAtoms.resize(local.num_local);
#endif

  if (comm.rank == 0) {
    fmt::printf("set global ids.\n");
    fflush(stdout);
  }

#ifdef LSMS_DEBUG
  MPI_Barrier(comm.comm);
  // write the distribution of atoms
  std::vector<int> atomsPerNode(comm.size);
#endif

  // set up exchange correlation functionals
  if (lsms.xcFunctional[0] == 1) {
    lsms.exch_corr =
        std::make_shared<lsms::XCLibxc>(lsms.n_spin_pola, lsms.xcFunctional);
  }                                 // use libxc functional
  if (lsms.xcFunctional[0] == 2) {  // use new LSMS functional
    lsms.newFunctional.init(lsms.n_spin_pola, lsms.xcFunctional);
  }

  AngularMomentumIndices::init(2 * crystal.maxlmax);
  SphericalHarmonicsCoeficients::init(2 * crystal.maxlmax);

  GauntCoeficients::init(lsms);
  IFactors::init(lsms, crystal.maxlmax);

#if defined(ACCELERATOR_CUDA_C) || defined(ACCELERATOR_HIP)
  deviceConstants.allocate();
#endif

  double timeBuildLIZandCommList = MPI_Wtime();
  if (lsms.global.iprint >= 0) {
    fmt::printf(
        "building the LIZ and Communication lists [buildLIZandCommLists]\n");
    fflush(stdout);
  }
  buildLIZandCommLists(comm, lsms, crystal, local);
  timeBuildLIZandCommList = MPI_Wtime() - timeBuildLIZandCommList;
  if (lsms.global.iprint >= 0) {
    fmt::printf("time for buildLIZandCommLists [num_local=%d]: %lf sec\n",
                local.num_local, timeBuildLIZandCommList);
    fflush(stdout);
  }

#ifdef LSMS_DEBUG
  MPI_Barrier(comm.comm);
#endif

  // initialize the potential accelerators (GPU)
  // we need to know the max. size of the kkr matrix to invert:
  // lsms.n_spin_cant*local.maxNrmat() which is only available after building
  // the LIZ

  acceleratorInitialize(lsms.n_spin_cant * local.maxNrmat(),
                        lsms.global.GPUThreads);
  local.tmatStore.pinMemory();
#if defined(ACCELERATOR_CUBLAS) || defined(ACCELERATOR_LIBSCI) || \
    defined(ACCELERATOR_CUDA_C) || defined(ACCELERATOR_HIP)
  // deviceStorage = allocateDStore();
  deviceStorage = new DeviceStorage;
#endif

  for (int i = 0; i < local.num_local; i++)
    local.atom[i].pmat_m.resize(lsms.energyContour.groupSize());

  // set maximal number of radial grid points and core states if reading from
  // bigcell file
  local.setMaxPts(lsms.global.iprpts);
  local.setMaxCore(lsms.global.ipcore);

  if (lsms.global.iprint >= 0) printLSMSGlobals(stdout, lsms);
  if (lsms.global.iprint >= 0) printLSMSSystemParameters(stdout, lsms);
  if (lsms.global.iprint >= 1) printCrystalParameters(stdout, crystal);
  if (lsms.global.iprint >= 0) printAlloyParameters(stdout, alloyDesc);

  if (lsms.global.iprint >= 1) {
    printCommunicationInfo(stdout, comm);
  }
  fflush(stdout);

#ifdef LSMS_DEBUG
  if (lsms.global.iprint >= 0) {
    fmt::printf(
        "Entering the Voronoi construction BEFORE loading the potentials.\n");
    fflush(stdout);
  }
  MPI_Barrier(comm.comm);
#endif

  if (lsms.use_voronoi == 1) {
    setupVorpol(lsms, crystal, local);
  } else {
    lsms::setupCell(lsms, crystal, local);
    lsms::printCell(lsms, crystal, local, comm);
  }

#ifdef LSMS_DEBUG
  if (lsms.global.iprint >= 0) {
    fmt::printf("Entering the LOADING of the potentials.\n");
    fflush(stdout);
  }
  MPI_Barrier(comm.comm);
#endif

  double timeLoadPotential = MPI_Wtime();
  loadPotentials(comm, lsms, crystal, local);
  timeLoadPotential = MPI_Wtime() - timeLoadPotential;

  if (lsms.global.iprint >= 0) {
    fmt::printf("time loadPotential: %lf sec\n\n", timeLoadPotential);

    fprintf(stdout, "LIZ for atom 0 on this node\n");
    printLIZInfo(stdout, local.atom[0]);
    if (local.atom[0].forceZeroMoment) {
      fprintf(stdout, "\nMagnetic moment of atom 0 forced to be zero!\n\n");
    }
  }

  // Read evec file if in is define
  if (lsms.infoEvecFileIn[0] != 0) {
    readInfoEvec(comm, lsms, crystal, local, lsms.infoEvecFileIn);
    if (lsms.global.iprint >= 0) {
      fprintf(stdout, "Evec are read from: %s\n", lsms.infoEvecFileIn);
    }
  }

  if (!alloyDesc.empty()) {
    if (lsms.global.iprint >= 0) {
      fmt::printf("Entering the LOADING of the alloy banks.\n");
      fflush(stdout);
    }
    loadAlloyBank(comm, lsms, alloyDesc, alloyBank);
  }

#ifdef LSMS_DEBUG
  if (lsms.global.iprint >= 0) {
    fmt::printf(
        "Entering the Voronoi construction AFTER loading the potentials.\n");
    fflush(stdout);
  }
  MPI_Barrier(comm.comm);
#endif

  if (lsms.use_voronoi == 1) {
    setupVorpol(lsms, crystal, local);
  } else {
    lsms::setupCell(lsms, crystal, local);
  }

#ifdef LSMS_DEBUG
  MPI_Barrier(comm.comm);
#endif

  // Generate new grids after new rmt is defined
  for (int i = 0; i < local.num_local; i++) {
    if (local.atom[i].generateNewMesh)
      interpolatePotential(lsms, local.atom[i]);
  }

  double timeCalculateVolumes = MPI_Wtime();
  calculateVolumes(comm, lsms, crystal, local);
  timeCalculateVolumes = MPI_Wtime() - timeCalculateVolumes;
  if (lsms.global.iprint >= 0) {
    fmt::printf("time calculateVolumes: %lf sec\n\n", timeCalculateVolumes);
  }

  // initialize Mixing
  double timeSetupMixing = MPI_Wtime();

  bool activate_spin_mixing = false;

  // Charge density mixer
  std::unique_ptr<lsms::AbstractMixer> chdMixer;
  // Spin density mixer
  std::unique_ptr<lsms::AbstractMixer> spdMixer;
  // Initial charge and spin density mixer
  std::unique_ptr<lsms::AbstractMixer> initMixer;
  // Potential Mixer
  std::unique_ptr<lsms::AbstractMixer> potMixer;

  lsms::printMixingParameters(mix, comm, lsms);

  // Charge density vector
  lsms::ChargeMixingVector chdMixVector;
  // Spin density vector
  lsms::SpinMixingVector spdMixVector;

  if (mix.n_init_iterations > 0) {
    chdMixVector = lsms::ChargeMixingVector(lsms, local, true);
    chdMixer = lsms::AbstractMixer::generateMixer(mix.init_mixer_type, chdMixVector,
                                                  mix.initMixingParameter);
  }


  // Potential mixer
  lsms::PotentialMixingVector potMixVector(lsms, local);
  potMixer = lsms::AbstractMixer::generateMixer(mix.pot_mixer_type, potMixVector,
                                                mix.potMixingParameter);

  timeSetupMixing = MPI_Wtime() - timeSetupMixing;
  if (lsms.global.iprint >= 0) {
    fmt::printf("time setupMixing: %lf sec\n\n", timeSetupMixing);
  }

  double timeCalculateMadelungMatrix = MPI_Wtime();

  lsms::calculateMultiMadelungMatrices(lsms, crystal, local);

  if (lsms.global.debug_madelung) {
    lsms::printMultiMadelungMatrices(lsms, local, comm);
  }

  timeCalculateMadelungMatrix = MPI_Wtime() - timeCalculateMadelungMatrix;

  if (comm.rank == 0) {
    fmt::printf("time calculateMultiMadelungMatrices: %lf sec\n\n",
                timeCalculateMadelungMatrix);
  }

  if (lsms.global.iprint >= 1) {
    printLocalTypeInfo(stdout, local);
  }

  if (lsms.n_spin_cant > 1) {
    for (int i = 0; i < local.num_local; i++) local.atom[i].get_b_basis();
  } else {
    for (int i = 0; i < local.num_local; i++) local.atom[i].reset_b_basis();
  }

#ifdef USE_PAPI
#define NUM_PAPI_EVENTS 2
  int hw_counters = PAPI_num_counters();
  if (hw_counters > NUM_PAPI_EVENTS) hw_counters = NUM_PAPI_EVENTS;
  int papi_events[NUM_PAPI_EVENTS];
  char *papi_event_name[] = {"PAPI_TOT_INS", "PAPI_FP_OPS"};
  // get events from names:
  for (int i = 0; i < NUM_PAPI_EVENTS; i++) {
    if (PAPI_event_name_to_code(papi_event_name[i], &papi_events[i]) != PAPI_OK)
      if (hw_counters > i) hw_counters = i;
  }
  long long papi_values[NUM_PAPI_EVENTS + 4];
  if (hw_counters > NUM_PAPI_EVENTS) hw_counters = NUM_PAPI_EVENTS;
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

  if (lsms.global.iprint >= 0) {
    fmt::printf("Total number of iterations:%d\n", lsms.nscf);
    fflush(stdout);
  }

  // timing
  double timeScfLoop = MPI_Wtime();
  double timeCalcChemPot = 0.0;
  double timeCalcPotential = 0.0;
  double timeCalcCharge = 0.0;
  double timeCalcMixPotential = 0.0;
  double timeCalcMixCharge = 0.0;

  int iterationStart = 0;
  int potentialWriteCounter = 0;

  FILE *kFile = nullptr;
  if (comm.rank == 0) {
    iterationStart = readNextIterationNumber("k.out");
    kFile = fopen("k.out", "a");
  }

  std::unique_ptr<lsms::Potential> potential =
      std::make_unique<lsms::ASAPotential>();

  // Charges on sites
  std::vector<Real> qsub(crystal.num_types, 0.0);

  if (lsms.pot_in_type == -1) {

    potential->calculatePotential(comm, lsms, local, crystal, qsub);

    // Initialize potentials and charge densities
    lsms::copyChargesAndPotential(lsms, local);

    // Check charge density after mixing
    lsms::checkRadialChargeDensity(lsms, local);

    if (comm.rank == 0) {
      fmt::printf("Initial MTZ: %20.9f\n", lsms.vmt);
    }

  }

  int iteration = 0;

  lsms::MixingParameter spdMixingParameter;
  unsigned int spd_mixer_type;
  bool change_mixer = false;

#ifdef USE_NVTX
  nvtxRangePushA("SCFLoop");
#endif
  for (iteration = 0; iteration < lsms.nscf && !(converged && energyConverged);
       iteration++) {
    oldTotalEnergy = lsms.totalEnergy;

    std::fill(qsub.begin(), qsub.end(), 0.0);

    if (lsms.global.iprint >= -1 && comm.rank == 0)
      fmt::printf("SCF iteration %d:\n", iteration);

    // Recalculate core states from `vr`
    lsms::calculateCoreStates(comm, lsms, local);

    // Calculate band energy from `vr`
    energyContourIntegration(comm, lsms, local);

    // Calculate chemical potential
    double dTimeCCP = MPI_Wtime();
    lsms::calculateChemPot(comm, lsms, local, eband);
    dTimeCCP = MPI_Wtime() - dTimeCCP;
    timeCalcChemPot += dTimeCCP;

    // Create radial charge density in `rhoNew`
    double dTimeCHD = MPI_Wtime();
    lsms::calculateRadialChargeDensity(lsms, local);

    // Calculate charge density rms `New`
    lsms::calculateQRMS(lsms, local);

    // Create charges
    lsms::calculateCharge(lsms, local, qsub);
    globalSum(comm, qsub.data(), crystal.num_types);

    dTimeCHD = MPI_Wtime() - dTimeCHD;
    timeCalcCharge += dTimeCHD;

    // Mix charge density
    double dTimeCM = MPI_Wtime();

    // Check charge density after mixing
    lsms::checkRadialChargeDensity(lsms, local);

    // Change representation
    if (lsms.n_spin_pola == 2) {

      for (auto i = 0; i < local.atom.size(); i++) {
        for (auto j = 0; j < local.atom[i].jmt; j++) {

          local.atom[i].rhoNew(j, 0) = local.atom[i].rhoNew(j, 0) + local.atom[i].rhoNew(j, 1);
          local.atom[i].rhoNew(j, 1) = local.atom[i].rhoNew(j, 0) - 2.0 * local.atom[i].rhoNew(j, 1);

          local.atom[i].rhotot(j, 0) = local.atom[i].rhotot(j, 0) + local.atom[i].rhotot(j, 1);
          local.atom[i].rhotot(j, 1) = local.atom[i].rhotot(j, 0) - 2.0 * local.atom[i].rhotot(j, 1);

        }
      }

    }

    // Init charge density or initial spin density mixer
    if (mix.n_init_iterations == iteration) {

      if (mix.n_init_spin_iterations > 0) {
        spdMixingParameter = mix.initSpdMixingParameter;
        spd_mixer_type = mix.init_spd_mixer_type;
        if (comm.rank == 0) {
          fmt::printf(" === Change Mixer (Initial) ===\n");
        }
      } else {
        spdMixingParameter = mix.spdMixingParameter;
        spd_mixer_type = mix.spd_mixer_type;
        if (comm.rank == 0) {
          fmt::printf(" === Change Mixer ===\n");
        }
      }

      change_mixer = true;

    }

    if (mix.n_init_spin_iterations > 0) {
      if (mix.n_init_spin_iterations == (iteration - mix.n_init_iterations)) {
        if (comm.rank == 0) {
          fmt::printf(" === Change Mixer (After) ===\n");
        }
        spdMixingParameter = mix.spdMixingParameter;
        spd_mixer_type = mix.spd_mixer_type;
        change_mixer = true;
      }
    }

    if (change_mixer) {

      if (spd_mixer_type == lsms::MixerType::NO_MIXER | lsms.n_spin_pola == 1) {

        chdMixVector = lsms::ChargeMixingVector(lsms, local, true);
        chdMixer = lsms::AbstractMixer::generateMixer(mix.chd_mixer_type, chdMixVector,
                                                      mix.chdMixingParameter);

      } else {

        chdMixVector = lsms::ChargeMixingVector(lsms, local, false);
        spdMixVector = lsms::SpinMixingVector(lsms, local);

        chdMixer = lsms::AbstractMixer::generateMixer(mix.chd_mixer_type, chdMixVector,
                                                      mix.chdMixingParameter);

        spdMixer = lsms::AbstractMixer::generateMixer(spd_mixer_type, spdMixVector,
                                                      spdMixingParameter);

        activate_spin_mixing = true;
      }

      change_mixer = false;

    }

    chdMixVector.copyToVector(comm, lsms, local);

    if (activate_spin_mixing)
      spdMixVector.copyToVector(comm, lsms, local);

    chdMixer->mix(comm, chdMixVector);

    if (activate_spin_mixing)
      spdMixer->mix(comm, spdMixVector);

    chdMixVector.copyFromVector(comm, lsms, local);

    if (activate_spin_mixing)
      spdMixVector.copyFromVector(comm, lsms, local);

    if (lsms.n_spin_pola == 2) {

      for (auto i = 0; i < local.atom.size(); i++) {
        for (auto j = 0; j < local.atom[i].jmt; j++) {

          local.atom[i].rhoNew(j, 0) = 0.5 * (local.atom[i].rhoNew(j, 0) + local.atom[i].rhoNew(j, 1));
          local.atom[i].rhoNew(j, 1) = local.atom[i].rhoNew(j, 0) - local.atom[i].rhoNew(j, 1);

          local.atom[i].rhotot(j, 0) = 0.5 * (local.atom[i].rhotot(j, 0) + local.atom[i].rhotot(j, 1));
          local.atom[i].rhotot(j, 1) = local.atom[i].rhotot(j, 0) - local.atom[i].rhotot(j, 1);

        }
      }

    }

    // Check charge density after mixing
    lsms::checkRadialChargeDensity(lsms, local);

    dTimeCM = MPI_Wtime() - dTimeCM;
    timeCalcMixCharge += dTimeCM;

    // Create new potential `vrNew` from charge `rhoNew`
    double dTimeP = MPI_Wtime();
    potential->calculatePotential(comm, lsms, local, crystal, qsub);

    // Calculate charge density rms
    lsms::calculateVRMS(lsms, local);

    // Calculate total energy from `rhoNew` and `vr`
    calculateTotalEnergy(comm, lsms, local, crystal);

    dTimeP = MPI_Wtime() - dTimeP;
    timeCalcPotential += dTimeP;

    // Mix charge density
    double dTimePM = MPI_Wtime();

    potMixVector.copyToVector(comm, lsms, local);
    potMixer->mix(comm, potMixVector);
    potMixVector.copyFromVector(comm, lsms, local);

    dTimePM = MPI_Wtime() - dTimePM;
    timeCalcMixPotential += dTimePM;

    // RMS of charge and potential
    Real qrms = 0.0;
    Real vrms = 0.0;

    for (int i = 0; i < local.num_local; i++) {
      if (lsms.n_spin_pola == 2) {
        qrms = std::max(qrms,
                        0.5 * (local.atom[i].qrms[0] + local.atom[i].qrms[1]));
        vrms = std::max(vrms,
                        0.5 * (local.atom[i].vrms[0] + local.atom[i].vrms[1]));
      } else {
        qrms = std::max(qrms, local.atom[i].qrms[0]);
        vrms = std::max(vrms, local.atom[i].vrms[0]);
      }
    }

    globalMax(comm, qrms);
    globalMax(comm, vrms);

    // check for convergence
    converged = qrms < lsms.rmsTolerance;

    if (lsms.energyTolerance > 0)
      energyConverged = std::abs((lsms.totalEnergy - oldTotalEnergy) /
          lsms.totalEnergy) < lsms.energyTolerance;
    else
      energyConverged = true;

    Real mag = 0.0;
    if (lsms.n_spin_pola == 2) {

      for (int i = 0; i < local.num_local; i++) {
        mag += local.atom[i].xvalwsNew[0] - local.atom[i].xvalwsNew[1];
      }

      globalSum(comm, mag);
      mag /= (Real) lsms.num_atoms;

    }

    if (comm.rank == 0) {
      int gap_size = 12;
      int size = -1;
      size = std::max(size, num_digits(static_cast<int>(eband)));
      size = std::max(size, num_digits(static_cast<int>(lsms.totalEnergy)));
      size = std::max(size, num_digits(static_cast<int>(lsms.chempot)));
      size = std::max(size, num_digits(static_cast<int>(lsms.vmt)));

      gap_size -= size;
      gap_size = std::max(gap_size, 2);
      size += 10;

      fmt::printf("MTZ          = %*.9f Ry\n", size, lsms.vmt);
      fmt::printf("Band Energy  = %*.9f Ry %*s Fermi Energy = %15.12f Ry\n",
                  size, eband, gap_size, "", lsms.chempot);
      fmt::printf("Total Energy = %*.9f Ry\n", size, lsms.totalEnergy);
      fmt::printf("QRMS         = %*.5e\n", size, qrms);
      fmt::printf("VRMS         = %*.5e\n", size, vrms);
      fmt::printf("atol. Energy = %*.5e\n", size,
                  std::abs((lsms.totalEnergy - oldTotalEnergy)));

      if (lsms.n_spin_pola == 2) {
        fmt::printf("ave. mag.    = %*.9f\n", size, mag);
      }

    }

    // Update potentials and charge density
    lsms::updateChargePotential(lsms, local);

    if (kFile != nullptr) {
      fmt::print(kFile, "{:4d} {:22.12f} {:15.12f} {:12.6e} {:12.6e}\n",
                 iterationStart + iteration, lsms.totalEnergy, lsms.chempot,
                 qrms, vrms);
      fflush(kFile);
    }

    // Periodically write the new potential for scf calculations
    potentialWriteCounter++;
    if ((lsms.pot_out_type >= 0 && potentialWriteCounter >= lsms.writeSteps) ||
        converged) {
      if (comm.rank == 0)
        std::cout << "Writing new potentials and restart file.\n";
      writePotentials(comm, lsms, crystal, local);
      potentialWriteCounter = 0;
      if (comm.rank == 0) {
        writeRestart("i_lsms_restart.lua", lsms, crystal, mix, potentialShifter,
                     alloyDesc);
      }
    }

  }  // end of iteration

#ifdef USE_NVTX
  nvtxRangePop();
#endif
  timeScfLoop = MPI_Wtime() - timeScfLoop;

// Recalculate core states
  lsms::calculateCoreStates(comm, lsms, local
  );

  writeInfoEvec(comm, lsms, crystal, local, eband, lsms
      .infoEvecFileOut);
  if (lsms.localAtomDataFile[0] != 0)
    writeLocalAtomData(comm, lsms, crystal, local, eband,
                       lsms
                           .localAtomDataFile);

  if (kFile != nullptr)
    fclose(kFile);

/**
 * Total energy calculation of all contributions
 */

  lsms::DFTEnergy dft_energy;
  calculateTotalEnergy(comm, lsms, local, crystal, dft_energy
  );

  if (comm.rank == 0) {
    lsms::print_dft_energy(dft_energy);
  }

// -----------------------------------------------------------------------------

#ifdef USE_PAPI
  PAPI_stop_counters(papi_values, hw_counters);
  papi_values[hw_counters] = PAPI_get_real_cyc() - papi_real_cyc_0;
  papi_values[hw_counters + 1] = PAPI_get_real_usec() - papi_real_usec_0;
  papi_values[hw_counters + 2] = PAPI_get_virt_cyc() - papi_virt_cyc_0;
  papi_values[hw_counters + 3] = PAPI_get_virt_usec() - papi_virt_usec_0;
  long long accumulated_counters[NUM_PAPI_EVENTS + 4];
  MPI_Reduce(papi_values, accumulated_counters, hw_counters + 4, MPI_LONG,
             MPI_SUM, 0, MPI_COMM_WORLD);
  if (comm.rank == 0) {
    for (int i = 0; i < hw_counters; i++)
      std::cout << "Accumulated: " << (papi_event_name[i]) << " = "
                << (accumulated_counters[i]) << "\n";
    std::cout << "PAPI accumulated real cycles : "
              << (accumulated_counters[hw_counters]) << "\n";
    std::cout << "PAPI accumulated user cycles : "
              << (accumulated_counters[hw_counters + 2]) << "\n";
    double gflops_papi = ((double)accumulated_counters[1]) /
                         (1000.0 * (double)papi_values[hw_counters + 1]);
    double gflops_hw_double = ((double)accumulated_counters[2]) /
                              (1000.0 * (double)papi_values[hw_counters + 1]);
    double gflops_hw_single = ((double)accumulated_counters[3]) /
                              (1000.0 * (double)papi_values[hw_counters + 1]);
    double gips = ((double)accumulated_counters[0]) /
                  (1000.0 * (double)papi_values[hw_counters + 1]);
    std::cout << "PAPI_FP_OPS real GFLOP/s : " << (gflops_papi) << "\n";
    std::cout << "PAPI hw double real GFLOP/s : " << (gflops_hw_double) << "\n";
    std::cout << "PAPI hw single real GFLOP/s : " << (gflops_hw_single) << "\n";
    std::cout << "PAPI real GINST/s : " << (gips) << "\n";
    std::cout << "Time (s) : " << (double)papi_values[hw_counters + 1] << "\n";
  }
#endif

  if (lsms.pot_out_type >= 0) {
    if (comm.rank == 0) std::cout << "Writing new potentials.\n";
    writePotentials(comm, lsms, crystal, local
    );
    if (comm.rank == 0) {
      std::cout << "Writing restart file.\n";
      writeRestart("i_lsms_restart.lua", lsms, crystal, mix, potentialShifter,
                   alloyDesc);
    }
  }

  double fomScale = calculateFomScaleDouble(comm, local);

  auto lsmsEndTime = std::chrono::steady_clock::now();
  std::chrono::duration<double> lsmsRuntime = lsmsEndTime - lsmsStartTime;
  std::chrono::duration<double> lsmsInitTime = lsmsEndInitTime - lsmsStartTime;

  if (comm.rank == 0) {
    int size = -1;
    size = std::max(size, num_digits(static_cast<int>(eband)));
    size = std::max(size, num_digits(static_cast<int>(lsms.chempot)));
    size = std::max(size, num_digits(static_cast<int>(lsms.totalEnergy)));
    size += 17;

    fmt::printf("Band Energy  = %*.15f Ry\n", size, eband);
    fmt::printf("Fermi Energy = %*.15f Ry\n", size, lsms.chempot);
    fmt::printf("Total Energy = %*.15f Ry\n", size, lsms.totalEnergy);

    fmt::printf("\nTimings:\n========\n");

    fmt::printf("LSMS Runtime               = %lf s\n", lsmsRuntime.
        count()
    );
    fmt::printf("LSMS Initialization Time   = %lf s\n", lsmsInitTime.
        count()
    );
    fmt::printf("timeScfLoop                = %lf s\n", timeScfLoop);
    fmt::printf("     number of iteration:%d\n", iteration);
    fmt::printf("timeScfLoop/iteration      = %lf s\n",
                timeScfLoop / (double) iteration);
    fmt::printf("timeCalcChemPot/iteration  = %lf s\n",
                timeCalcChemPot / (double) iteration);
    fmt::printf("timeCalcPotenial/iteration = %lf s\n",
                timeCalcPotential / (double) iteration);
    fmt::printf("timeCalcCharge/iteration   = %lf s\n",
                timeCalcCharge / (double) iteration);
    fmt::printf("timeCalcMixPotal/iteration = %lf s\n",
                timeCalcMixPotential / (double) iteration);
    fmt::printf("timeCalcMixChd/iteration   = %lf s\n",
                timeCalcMixCharge / (double) iteration);
    fmt::printf("timeBuildLIZandCommList    = %lf s\n",
                timeBuildLIZandCommList);

    long energyContourPoints = 1;
    if (lsms.energyContour.grid == 2) {
      energyContourPoints = lsms.energyContour.npts + 1;
    }
    fmt::printf("FOM Scale                  = %lf\n", (double) fomScale);
    fmt::printf("Energy Contour Points      = %ld\n", energyContourPoints);
    fmt::printf("FOM / energyContourPoint   = %lg/s\n",
                fomScale * (double) iteration / timeScfLoop);
    fmt::printf("FOM                        = %lg/s\n",
                (double) energyContourPoints * (double) fomScale *
                    (double) iteration / timeScfLoop);
  }

  local.tmatStore.
      unpinMemory();

#if defined(ACCELERATOR_CUBLAS) || defined(ACCELERATOR_LIBSCI) || \
    defined(ACCELERATOR_CUDA_C) || defined(ACCELERATOR_HIP)
  // freeDStore(deviceStorage);
  delete deviceStorage;
#endif

  acceleratorFinalize();

#ifdef USE_GPTL
  GPTLpr(comm.rank);
#endif

  H5close();
  finalizeCommunication();
  lua_close(L);
  return EXIT_SUCCESS;
}
