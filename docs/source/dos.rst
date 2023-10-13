******************
Density of States
******************

To get DOS, first a standard self-consistent calculation (as described in the previous sections) must be performed. Once that is done, copy the restart file to a new input file and the converged potential as the new starting potential

.. parsed-literal::
   cp i_lsms.restart i_lsms.dos
   cp w_<systemid> v_<systemid>

Edit the i_lsms.dos for dos calculations. First set the mode for DOS

.. parsed-literal::
   lsmsMode="dos"

Change the energy contour settings so that the points are parallel to the real energy axis

.. parsed-literal::
   energyContour.grid=3
   energyContour.npts= <number of energy points>
   energyContour.ebot= <starting energy value>
   energyContour.etop= <final energy value>
   energyContour.eibot= <imaginary part>

It is important to specify ``energyContour.eibot`` as we cannot traverse the real energy line due to singularities. Hence adding a small imaginary part to the energy is necessary. Generallyvalues between 0.001 and 0.005 Ryd are good choices.

lsms will generate an output file called ``dos.out`` with three columns

1. Real part of the energy.
2. Imaginary part of the energy.
3. Total Density of States
