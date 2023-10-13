******************************************
Green's Function at Matsubara Frequencies
******************************************

lsms can write out the Green's function at the Matsubara frequencies. To do this, set the mode and the required temperature

.. parsed-literal::
   lsmsMode = "gf_out"
   temperature = <Temperature in Kelvin>

The relevant energy points are given by

.. math::
   E_n = E_F + \\pi ik_B T(2n + 1)

To generate these energy points, we need to set the appropriate contour in the input file

.. parsed-literal::
   energyContour.grid=4
   energyContour.npts= <number of points>

lsms will generate output files of the form ``greens_functions_<atom index>.out``, containing the local Green's function in blocks for each energy. There are five columns in each energy block

1. Spin index
2. L = (l,m) index
3. L' = (l',m') index
4. Real part of G(L,L')
5. Imaginary part of G(L,L') 
