**********
Input File
**********

The input file is a script in the LUA language (lua.org)
lsms will read the input values from LUA variables or use default values when appropriate. 
Lines starting with -- are comments. 

Let's go through the different variables that need to be specified in the input file for the example system FeCo

.. parsed-literal::
   systemid="FeCo"
   system_title="Iron-Cobalt test for LSMS 3"

``systemid`` is used to construct potential filenames.
``system_title`` is an arbitrary text to describe the system.

.. parsed-literal::
   pot_in_type=1
   pot_out_type=0

``pot_in_type`` and ``pot_out_type`` are the file formats for the potential files. There are three options available

1. ``...type=0`` - HDF5 binary file. Only one file is needed as input and only one file will be produced as output for the whole system.
2. ``...type=1`` - Text file. One file per atom is needed/produced. For eg, if there are 100 atoms in the system and type=1 is chosen, then 100 starting potential files are needed and 100 converged potential files will be produced.
3. ``...type=-1`` - No starting potential needed and no output potential file will be produced. The code will generate a starting potential.

.. parsed-literal::
   num_atoms=2
   nspin=3
   nscf=10

``num_atoms`` represents the number of atoms in the simulation cell.
``nspin`` determines how magnetism is treated. There are three options

1. ``nspin=1`` : no spin polarization
2. ``nspin=2`` : collinear spins
3. ``nspin=3`` : non-collinear magnetism

``nscf`` is the maximum number of self-consistent iterations.

.. parsed-literal::
   mixing = { {quantity = "potential", algorithm = "broyden", mixing_parameter = 0.05} }
   numberOfMixQuantities = 0
   for k,v in pairs(mixing) do
     numberOfMixQuantities = numberOfMixQuantities + 1
   end

This sets the mixing of quantities during the self-consistent iterations, generally this does not need to be modified. The last four lines are boilerplate to count the number of quantities to be mixed. 

.. parsed-literal::
   energyContour = {npts=31, grid=2, ebot=-0.3, etop=0.0, eitop=0.825, eibot=0.0025}

This specifies the contour for integration of the Green's function. ``npts`` is the number of grid points along the contour. ``ebot`` is the bottom of the contour, ``etop`` is the top of the contour which is always set to the Fermi energy (``etop = 0``). ``eitop`` and ``eibot`` are the imaginary parts of the energy at the top and end of the contour. 

.. parsed-literal::
   site_default = {lmax=3, rLIZ=13.5, rsteps={96.9, 97.9, 98.9, 99.9}}

``lmax`` represents the angular momentum cutoff and ``rLIZ`` is the radius of the local interaction zone. Refer to the theory behind LSMS to know more about the local interaction zone. 


