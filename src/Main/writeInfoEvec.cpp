/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */

#include "writeInfoEvec.hpp"

// write out the magnetism and constraint info for each site
static void writeSingleEvec(FILE *f,int z, int i, Real posX, Real posY, Real posZ, AtomData &atom)
{
// Z global_id x y z  qtotws  mtotws  evec_x evec_y evec_z  e_mix  B_x B_y B_z  vSpinShift
  fprintf(f,"%3d %8d  %21.15lf  %21.15lf  %21.15lf  %12.6lf %12.6lf  %21.15lf  %21.15lf  %21.15lf  %6.2lf  %21.15lf  %21.15lf  %21.15lf  %8.4lf  %21.15lf  %21.15lf  %21.15lf\n",
          z,i, posX, posY, posZ,
          atom.qtotws, atom.mtotws,
          atom.evec[0], atom.evec[1], atom.evec[2],
          -1.0,
          atom.b_con[0], atom.b_con[1], atom.b_con[2],
          atom.vSpinShift,
          atom.evecOut[0], atom.evecOut[1], atom.evecOut[2]);
}

// write out the magnetism and constraint info for each site
static void writeSingleLocalAtomData(FILE *f,int z, int i, Real posX, Real posY, Real posZ, AtomData &atom)
{
  // Z global_id x y z  qtotws  mtotws  evec_x evec_y evec_z  B_x B_y B_z  vSpinShift localVolume localEnergy  mtotmt  mvalmt  mvalws  rmt jws
  fprintf(f,"%3d %8d  %21.15lf  %21.15lf  %21.15lf  %12.6lf %12.6lf  %21.15lf  %21.15lf  %21.15lf  %21.15lf  %21.15lf  %21.15lf  %8.4lf  %.15lf  %.15lf %21.15lf %21.15lf %21.15lf %10.6lf %10.6lf\n",
          z,i, posX, posY, posZ,
          atom.qtotws, atom.mtotws,
          atom.evec[0], atom.evec[1], atom.evec[2],
          atom.b_con[0], atom.b_con[1], atom.b_con[2],
          atom.vSpinShift,
          atom.omegaWS, atom.localEnergy+atom.localMadelungEnergy,
          atom.mtotmt, atom.mvalmt, atom.mvalws,
          atom.r_mesh[atom.jmt], atom.r_mesh[atom.jws]
          );
}

static void readSingleEvec(FILE *f,int &z, int &i, Real &posX, Real &posY, Real &posZ, AtomData &atom)
{
  Real tmp1, tmp2;
// Z global_id x y z  qtotws  mtotws  evec_x evec_y evec_z  e_mix  B_x B_y B_z  vSpinShift
  auto retval = fscanf(f,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n"
         ,&z,&i, &posX, &posY, &posZ, &tmp1, &atom.mtotws,
         &atom.evec[0], &atom.evec[1], &atom.evec[2],
         &tmp2,
         &atom.b_con[0], &atom.b_con[1], &atom.b_con[2],
         &atom.vSpinShift,
         &atom.evecOut[0], &atom.evecOut[1], &atom.evecOut[2]);
}

static void readSingleEvec(FILE *f, int &z, int &i,
                           Real *pos, Real &mtotws, Real *evec, Real *b_con, Real &vSpinShift,
                           Real *evecOut)
{
  Real tmp1, tmp2;
  // Z global_id x y z  qtotws  mtotws  evec_x evec_y evec_z  e_mix  B_x B_y B_z  vSpinShift
  auto retval = fscanf(f,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n"
         ,&z,&i, &pos[0], &pos[1], &pos[2], &tmp1, &mtotws,
         &evec[0], &evec[1], &evec[2],
         &tmp2,
         &b_con[0], &b_con[1], &b_con[2],
         &vSpinShift,
         &evecOut[0], &evecOut[1], &evecOut[2]);
}


int writeInfoEvec(LSMSCommunication &comm,LSMSSystemParameters &lsms, CrystalParameters &crystal, LocalTypeInfo &local, Real eband, const char *name)
{
  AtomData pot_data;

  if(comm.rank==0)
  {
    FILE *outf=fopen(name,"w");

    fprintf(outf,"%.15lf %.15lf %.15lf\n", lsms.totalEnergy, eband, lsms.chempot);
    
// loop over all atom types:
    for(int i=0; i<crystal.num_types; i++)
    {
      if(crystal.types[i].node==comm.rank)
      {
        writeSingleEvec(outf,crystal.types[i].Z,i,
                        crystal.position(0,crystal.types[i].first_instance), // posX
                        crystal.position(1,crystal.types[i].first_instance), // posY
                        crystal.position(2,crystal.types[i].first_instance), // posZ
                        local.atom[crystal.types[i].local_id]);
      } else {
        int local_id;
        communicateSingleAtomData(comm, crystal.types[i].node, comm.rank, local_id, pot_data,i);
        if(local_id!=crystal.types[i].local_id) printf("WARNING: local_id doesn't match in writePotentials!\n");
        writeSingleEvec(outf,crystal.types[i].Z,i,
                        crystal.position(0,crystal.types[i].first_instance), // posX
                        crystal.position(1,crystal.types[i].first_instance), // posY
                        crystal.position(2,crystal.types[i].first_instance), // posZ
                        pot_data);
      }
    }
    fclose(outf);
  } else { // comm.rank!=0
    for(int i=0; i<local.num_local; i++)
    {
      communicateSingleAtomData(comm, comm.rank, 0, i, local.atom[i],local.global_id[i]);
    }
  }

  return 0;
}

int readInfoEvec(LSMSCommunication &comm,LSMSSystemParameters &lsms, CrystalParameters &crystal, LocalTypeInfo &local, const char *name)
{
  AtomData pot_data;

  int Z,ii;
  Real pos[3], evec[3], b_con[3], evec_out[3];
  Real ef, etot, eband;
  Real mtotws, spin_shift;

  if(comm.rank==0)
  {
    FILE *inf=fopen(name,"r");
    auto retval = fscanf(inf,"%lf %lf %lf\n", &etot, &eband, &ef);
    
// loop over all atom types:
    for(int i=0; i<crystal.num_types; i++)
    {

      if(crystal.types[i].node==comm.rank)
      {
        // Atom types that are on the same process as the reading proce

        readSingleEvec(inf,Z,ii,
                       pos, mtotws, evec, b_con, spin_shift, evec_out);

        for (int j = 0; j < 3; j++) {
          local.atom[crystal.types[i].local_id].b_con[j] = b_con[j];
        }

// still need to set evec in crystal!
      } else {

        readSingleEvec(inf,Z,ii,
                       pos, mtotws, evec, b_con, spin_shift, evec_out);


        MPI_Send(b_con, 3, MPI_DOUBLE, crystal.types[i].node, ii, comm.comm);

      }
    }
    fclose(inf);
  } else { // comm.rank!=0

    for(int i=0; i<local.num_local; i++)
    {

      MPI_Status status;
      MPI_Recv(b_con, 3, MPI_DOUBLE, 0, local.global_id[i], comm.comm, &status);

      for (int j = 0; j < 3; j++) {
        local.atom[i].b_con[j] = b_con[j];
      }

    }
  }

  return 0;
}

 int writeLocalAtomData(LSMSCommunication &comm,LSMSSystemParameters &lsms, CrystalParameters &crystal, LocalTypeInfo &local, Real eband, const char *name)
{
  AtomData pot_data;

  if(comm.rank==0)
  {
    FILE *outf=fopen(name,"w");

    fprintf(outf,"# totalEnergy  bandEnergy  fermiEnergy electrostaticEnergy\n");
    fprintf(outf,"# Z global_id x y z  qtotws  mtotws  evec_x evec_y evec_z  e_mix  B_x B_y B_z  vSpinShift localVolume localEnergy  mtotmt  mvalmt  mvalws\n");
    
    fprintf(outf,"%.15lf %.15lf %.15lf %.15lf\n", lsms.totalEnergy, eband, lsms.chempot, lsms.u0);
    
// loop over all atom types:
    for(int i=0; i<crystal.num_types; i++)
    {
      if(crystal.types[i].node==comm.rank)
      {
        writeSingleLocalAtomData(outf,crystal.types[i].Z,i,
                        crystal.position(0,crystal.types[i].first_instance), // posX
                        crystal.position(1,crystal.types[i].first_instance), // posY
                        crystal.position(2,crystal.types[i].first_instance), // posZ
                        local.atom[crystal.types[i].local_id]);
      } else {
        int local_id;
        communicateSingleAtomData(comm, crystal.types[i].node, comm.rank, local_id, pot_data,i);
        if(local_id!=crystal.types[i].local_id) printf("WARNING: local_id doesn't match in writePotentials!\n");
        writeSingleLocalAtomData(outf,crystal.types[i].Z,i,
                        crystal.position(0,crystal.types[i].first_instance), // posX
                        crystal.position(1,crystal.types[i].first_instance), // posY
                        crystal.position(2,crystal.types[i].first_instance), // posZ
                        pot_data);
      }
    }
    fclose(outf);
  } else { // comm.rank!=0
    for(int i=0; i<local.num_local; i++)
    {
      communicateSingleAtomData(comm, comm.rank, 0, i, local.atom[i],local.global_id[i]);
    }
  }

  return 0;
}

