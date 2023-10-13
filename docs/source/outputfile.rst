*************
Output Files
*************

There are three output files/sources, the standard output, ``k.out`` and ``info_evec_out``

Standard Output
#################

The standard output (stdout) contains the full details of the calculation. If using a job scheduler like SLURM, this will be the slurm output file (``slurm-<job number>.out``).

The starting few lines of the stdout looks like the following

.. parsed-literal::
  LSMS_3: Program started
  Using 16000 MPI processes
  Using 7 OpenMP threads
  Found 1 HIP GPUs.
  Device 0:
  Reading input file 'i_lsms'
  ...

and the last few lines look like the following

.. parsed-literal::
   Timings:
   ========
   LSMS Runtime = 227.354840 sec
   LSMS Initialization Time = 47.622071 sec
   timeScfLoop[rank==0] = 169.298811 sec
     number of iteration:5
   timeScfLoop/iteration = 33.859762 sec
   timeCalcChemPot[rank==0]/iteration = 1.263892 sec
   timeCalcPotentialsAndMixing[rank==0]/iteration = 0.003911 sec
   timeBuildLIZandCommList[rank==0]: 5.613144 sec
   FOM Scale = 58892754944000000.000000
   Energy Contour Points = 31
   FOM / energyContourPoint = 1.73931e+15/sec
   FOM = 5.39187e+16/sec


Summarized Output (k.out)
#########################

The ``k.out`` file presents only certain data at every SCF iteration, presented as a list. It looks like the following

.. parsed-literal::
   0 -24348906.123388141394     0.660632     0.000000    0.7799487815
   1 -24348844.791788913310     0.660562     0.000000    0.0002161170
   2 -24348843.298677124083     0.660561     0.000000    0.0000038077
   3 -24348843.281400021166     0.660561     0.000000    0.0000000438
   4 -24348843.281203657389     0.660561     0.000000    0.0000000005

It has 5 columns.
1. Iteration number
2. Total energy in atomic units, or Rydbergs (1 Ry = 13.6 eV)
3. Fermi energy (in Ry)
4. Magnetic moment on atom 1 (in Bohr magnetons)
5. RMS error (representing convergence)

Atomic Data (info_evec_out)
############################

The ``info_evec_out`` file presents the charge and magnetic data for each atom. It looks like the following

.. parsed-literal::
   -24348843.281203657388687 63963.682064113578235 0.660561420034407
   13        0      0.000000000000000     35.117180099999999     17.558590049999999     13.333009     0.000000      
   0.000000000000000      0.000000000000000      1.000000000000000 ...
   13        1     87.792950249999990     11.705726700000000     58.528633499999998     13.325152     0.000000      
   0.000000000000000      0.000000000000000      1.000000000000000 ...

The first line contains the total energy, band energy and Fermi energy. Each subsequent line has 18 columns.
Column 1: Atomic number
Column 2: Site index (starting at 0)
Columns 3-5: Site coordinates (x,y,z)
Column 6: Charge at site
Column 7: Magnetic moment at site
Column 8-10: Direction of magnetic moment
Column 11-18: More information on magnetic calculation like torque, desired moment direction etc.
