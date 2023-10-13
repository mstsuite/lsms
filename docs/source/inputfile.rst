**********
Input File
**********

The input file is a script in the LUA language (lua.org)
lsms will read the input values from LUA variables or use default values when appropriate. 
Lines starting with -- are comments. 


Run Parameters
###############

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


Position Data
##################

In lsms, the position data should be provided in the input file. The lattice vectors can be defined the following way

.. parsed-literal::
   a = 5.218
   bravais = {}
   bravais[1]={a,0,0}
   bravais[2]={0,a,0}
   bravais[3]={0,0,a}

The units are Bohr radii, also known as atomic units (1 a.u. = 0.5291 :math:`\angstrom`).
Then to setup the cell some boilerplate code is needed

.. parsed-literal::
   site = {}
   for i=1,num_atoms do site[i]={} end

After this the atoms have to be defined.

.. parsed-literal::
   -- FIRST ATOM
   site[1].pos={0,0,0}
   site[1].evec={0,0,1}
   site[1].pot_in_idx=0
   site[1].atom="Fe"
   site[1].Z=26
   site[1].Zc=10
   site[1].Zs=8
   site[1].Zv=8
   
   -- SECOND ATOM
   site[2].pos={0.5*a,0.5*a,0.5*a}
   site[2].evec={0,0,1}
   site[2].pot_in_idx=1
   site[2].atom="Co"
   site[2].Z=27
   site[2].Zc=10
   site[2].Zs=8
   site[2].Zv=9

Here are the definitions for the different keywords

1. ``pos`` represents the position of the atom in atomic units
2. ``evec`` sets the direction of the spin quantization axis
3. ``pot_in_idx`` is the index of the potential input file. If ``pot_in_idx=0``, the input potential for that atom will be read from the file v_FeCo.0. 
4. ``atom`` is the name of the atomic species
5. ``Z`` is the atomic number
6.  ``Zc``, ``Zs`` and ``Zv`` are the number of core, semicore and valence electrons respectively. Note that these should sum to ``Z``.

Finally some additional boilerplate code is needed to copy values defined in ``site_default`` into the atomic sites that have not defined them.

.. parsed-literal::
   -- set site defaults
   for i=1,num_atoms do
    for k,v in pairs(site_default) do
     if(site[i][k]==nil) then site[i][k]=v end
    end
   end

Reading from Position File
###########################

Alternatively, it is possible to create a separate position file and have the input file read from it. This is a more convenient option when the number of atoms is very large. the position file can be formatted in any way the user desires, there is no fixed way. Appropriate parsing code should then be provided in the input file. Additionally, the position file can also be given any file name. For example, consider the following format for the positions.

.. parsed-literal::
   117.057267       0.0       0.0
   0.0       117.057267       0.0
   0.0       0.0       117.057267
   Al       0.0       35.1171801       17.55859005
   Al       87.79295024999999       11.7057267       58.5286335
   Al       14.632158375       90.719381925       67.307928525
   Al       11.7057267       23.4114534       76.08722355
   Al       32.190748424999995       67.307928525       20.485021725
   Al       58.5286335       11.7057267       87.79295024999999
   ...
   ...
   ...

The first three lines are the lattice vectors and subsequent lines are the atom species name and the position, expressed in Cartesian coordinates. Let's call this file ``position.dat``. To parse this, the following lines are added to the input file

.. parsed-literal::
   posfile = io.open("position.dat")
   bravais = {}
   for i = 1, 3 do
     l = next_line(posfile)
     x, y, z = l:match("([+-]?%d*%.%d*)%s+([+-]?%d*%.%d*)%s+([+-]?%d*%.%d*)")
     bravais[i] = {x, y, z}
   end

   for i = 1, num_atoms do
     l = next_line(posfile)
     at, x, y, z = l:match("(%a+)%s+([+-]?%d*%.%d*)%s+([+-]?%d*%.%d*)%s+([+-]?%d*%.%d*)")
     print(l)
     print(at, x, y, z)
     site[i].atom = at
     site[i].pos = {x, y, z}
   end

   posfile:close()

   -- set atom types
   for i = 1,num_atoms do
     atom_name = site[i].atom
     if(atom_name~=nil) then
       for k,v in pairs(atom_type[atom_name]) do
         site[i][k]=v
       end
     end
   end

   -- set site defaults
   for i =1,num_atoms do
     for k,v in pairs(site_default) do
       if(site[i][k]==nil) then site[i][k]=v end
     end
   end

This type of procedure can be carried out for any position format.

Restarting Calculations
########################

After a calculation is complete lsms will generate output potentials with filenames starting with ``w_`` and the restart input file ``i_lsms.restart``. To restart, copy the new potential to the old potential and run the new input file

.. parsed-literal::
   cp w_FeCo v_FeCo
   mpirun -np <number of MPI ranks> $LSMS_PATH/lsms i_lsms.restart
