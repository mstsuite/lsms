# Implementation Notes for DMFT

## What needs to be changed in LSMS

### Energy Grid Construction

### `tmatStore`



## What needs to be calculated

###Matsubara Frequencies

For DMFT on Matsubara frequencies we will need to provide the local Green's function at the Matsubara frequencies $i\omega_n$.

$$\omega_n = \pi(2n+1)/\beta$$

with $\beta$ the inverse temperature

$$\beta = (k_B T)^{-1}$$

The Boltzmann constant is
$$
\begin{align}
k_B &= 1.380\,648\,52\times 10^{-23}\mathrm{J/K}\\
&= 8.617\,333\,262\,145\times 10^{-5}\mathrm{eV/K}\\
&= 6.333\,623\,179\,977\times 10^{-6}\mathrm{Ry/K}
\end{align}
$$


### Single Scattering Solutions

Calculate the regular and irregular Single Scatterer solutions $Z_i(r;\epsilon)$ and $J_i(r,\epsilon)$ and scattering $t_i(\epsilon)$ matrices for the Matsubara frequencies.

### Multiple Scattering Solution

Calculate the local scattering path matrix $\tau_{ii}(\epsilon)$

### Local Green's Function

Calculate the local Green's Function
$$
\begin{align}
G_{ii}(\epsilon) &=\int_{\Omega_i}G(r,r;\epsilon)d^3r\\
&= \tau_{ii}\int_{\Omega_i}Z(r;\epsilon)Z^+(r;\epsilon)d^3r - \int_{\Omega_i}Z(r;\epsilon)J^+(r;\epsilon)d^3r
\end{align}
$$
`localG(L,L',e)`

## New routines

- [ ] `void generateMatsubaraFrequencies(realPart, temparature, number, energies)`
  * `Real realPart` the real part of the energies generated
  * `Real temperature` the temperature in Kelvin
  * `int number` the number of energies to be generated
  * `std::vector<Complex> &energies` the vector of the `number` generated energies: `energies[n]`=`realPart`+$i\omega_n$
- [ ] `void calculateSingleScattererSolutions()`
- [ ] `void calculateLocalGreensFunctions()`
- [ ] 