/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */

#include <cstdio>
#include "AtomData.hpp"

int printAtomPotential(FILE *of, AtomData &atom)
{
  fprintf(of,"# Interstitial xc potentials and energies:\n");
  fprintf(of,"# exchangeCorrelationV[] = %lg %lg\n",atom.exchangeCorrelationV[0],atom.exchangeCorrelationV[1]);
  fprintf(of,"# exchangeCorrelationE = %lg\n", atom.exchangeCorrelationE);
  fprintf(of,"# vrms[] =  %lg %lg\n", atom.vrms[0], atom.vrms[1]);


  for(int ir=0; ir<atom.r_mesh.size(); ir++)
    fprintf(of,"%d %.10lg   %.10lg %.10lg   %.10lg %.10lg   %.10lg %.10lg  %.10lg %.10lg   %.10lg %.10lg  %.10lg %.10lg\n",
            ir,atom.r_mesh[ir],
            atom.rhoNew(ir,0), atom.rhoNew(ir,1),
            atom.rhotot(ir,0), atom.rhotot(ir,1), 
            atom.exchangeCorrelationPotential(ir,0), atom.exchangeCorrelationPotential(ir,1),
            atom.exchangeCorrelationEnergy(ir,0), atom.exchangeCorrelationEnergy(ir,1),
            atom.vr(ir,0), atom.vr(ir,1),
            atom.vrNew(ir,0), atom.vrNew(ir,1));

  // printf("Wrote vr_pot.out.\n");

  return 0;
}

