************************
Electrical Resistivity
************************

lsms can also calculate electrical resistivity from the Kubo-Greenwood equation. To do this, first perform a standard self-consistent calculation. Then copy the converged potential to the starting potential and copy the restart file

.. parsed-literal::
   cp w_<systemid> v_<systemid>
   cp i_lsms.restart i_lsms.dos

Change the mode to conductivity

.. parsed-literal::
   lsmsMode="kubo"

Then run lsms

.. parsed-literal::
   mpirun -np <number of MPI ranks> $LSMS_PATH/lsms i_lsms.kubo

The stdout will print the electrical resistivity in units of :math:`\mu\Omega`-cm

.. parsed-literal::
   TOTAL RESISTIVITY (in microOhm-cm)
   412.992   -1.32721   -1.86057
   0.540036   408.027   -0.912385
   -2.64419   -0.598773   416.019

For a spin-polarized system, the spin up and down resistivities are calculated separately and the resistors are added in parallel.

.. parsed-literal::
   RESISTIVITY SPIN 1 (in microOhm-cm)
   289.029  -0.206822  0.193379
   -0.222593  289.206  0.0590624
   0.0231963  -0.48124  289.46

   RESISTIVITY SPIN 2 (in microOhm-cm)
   262.711  0.162251  -0.491042
   -0.173043  261.377  0.476058
   -0.993332  1.08469  262.733

   TOTAL RESISTIVITY (in microOhm-cm)
   137.621   -0.00159926   -0.0911053
   -0.0977167   137.294   0.144432
   -0.26754   0.190146   137.724

Currently the resistivity implementation is non-relativistic. A relativistic implementation will be added in the future. 
