************
Running LSMS
************

Two types of files are necessary to run LSMS

1. **Input file** - This will contain all the calculation parameters, atom and crystal configuration. The default name for this file is i_lsms, although any file name can be chosen. 
2. **Starting potentials** - Need to be specified for each atom. 

If the binary has been compiled in $LSMS_PATH, the code can be run using

.. parsed-literal::
   mpirun -np <number of MPI ranks> $LSMS_PATH/lsms <name_of_input_file>

If the input filename is kept as i_lsms (the default), you can also run

.. parsed-literal::
   mpirun -np <number of MPI ranks> $LSMS_PATH/lsms 

Examples for runnning LSMS are present in lsms/Test. A few examples are also present in the documentation.
