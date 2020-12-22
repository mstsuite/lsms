/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#include "Main/SystemParameters.hpp"
#include "Communication/LSMSCommunication.hpp"
#include "CoreStates.hpp"

const bool useNewGetCoreStates =  true;

void calculateCoreStates(LSMSCommunication &comm, LSMSSystemParameters &lsms, LocalTypeInfo &local)
{
  for(int i=0; i<local.num_local; i++)
  {
    // dprintf("-- LSMS %d: Calculating core levels for local atom %d.\n", comm.rank, i);
    if(lsms.global.iprint>=0) printf("\ncalculateCoreState %d.%d\n",comm.rank,i);
    if(useNewGetCoreStates)
    {
      getCoreStates(lsms, local.atom[i]);
    } else {
      int local_iprpts=local.atom[i].vr.l_dim();
      int local_ipcore=local.atom[i].ec.l_dim();

      getcor_(&lsms.n_spin_pola,&lsms.mtasa,
	      &local.atom[i].jmt,&local.atom[i].jws,&local.atom[i].r_mesh[0],&local.atom[i].h,&local.atom[i].xstart,
	      &local.atom[i].vr(0,0),
	      &local.atom[i].numc,&local.atom[i].nc(0,0),&local.atom[i].lc(0,0),&local.atom[i].kc(0,0),&local.atom[i].ec(0,0),
	      &local.atom[i].ztotss,&local.atom[i].zsemss,&local.atom[i].zcorss,
	      &local.atom[i].ecorv[0],&local.atom[i].esemv[0],&local.atom[i].corden(0,0),&local.atom[i].semcor(0,0),
	      &lsms.nrelc,
	      &local.atom[i].qcpsc_mt,&local.atom[i].qcpsc_ws,&local.atom[i].mcpsc_mt,&local.atom[i].mcpsc_ws,
	      &local_iprpts,&local_ipcore,
	      &lsms.global.iprint,lsms.global.istop,32);
      local.atom[i].movedToValence[0] = local.atom[i].movedToValence[1] = 0;
      for(int j=0; j<local.atom[i].numc; j++)
	{
	  local.atom[i].coreStateType(j, 0) = 'C';
	  local.atom[i].coreStateType(j, 1) = 'C';
	}
    }
  }

// calculate the global maximum of ec:
  Real etopcor=-10.0e+20;
  for(int i=0; i<local.num_local; i++)
  {
    if(local.atom[i].numc>0)
    {
      for(int ic=0; ic<local.atom[i].numc; ic++)
        etopcor=std::max(local.atom[i].ec(ic,0),etopcor);
      if(local.atom[i].nspin>1) for(int ic=0; ic<local.atom[i].numc; ic++)
        etopcor=std::max(local.atom[i].ec(ic,1),etopcor);
    }
  }
  globalMax(comm,etopcor);

  lsms.largestCorestate=etopcor;
  if(lsms.global.iprint>=0)
    printf("Maximal Core State = %gRy\n",lsms.largestCorestate);
/*
      if(etopcor+0.1d0 .gt. ebot) then
         write(6,'('' GETCOR: etopcor+0.1 .gt. ebot'',2d12.5)')
     >         etopcor,ebot
         call kill(0,9)
         call fstop(sname)
      endif
*/
}

void calculateCoreStates(LSMSCommunication &comm, LSMSSystemParameters &lsms, AlloyAtomBank &alloyBank) {

  // for Wang-Landau for metallic alloys
  // overloaded existing routine for alloy potentials

  for(int i = 0; i < alloyBank.size(); i++)
  for(int j = 0; j < alloyBank[i].size(); j++) {
    AtomData &atom = alloyBank[i][j];

    int local_iprpts = atom.vr.l_dim();
    int local_ipcore = atom.ec.l_dim();
    if(lsms.global.iprint>=0) printf("\nalloy bank : calculateCoreState %d.%d\n",comm.rank,i);
    getcor_(&lsms.n_spin_pola,&lsms.mtasa,
            &atom.jmt,&atom.jws,&atom.r_mesh[0],&atom.h,&atom.xstart,
            &atom.vr(0,0),
            &atom.numc,&atom.nc(0,0),&atom.lc(0,0),&atom.kc(0,0),&atom.ec(0,0),
            &atom.ztotss,&atom.zsemss,&atom.zcorss,
            &atom.ecorv[0],&atom.esemv[0],&atom.corden(0,0),&atom.semcor(0,0),
            &lsms.nrelc,
            &atom.qcpsc_mt,&atom.qcpsc_ws,&atom.mcpsc_mt,&atom.mcpsc_ws,
            &local_iprpts,&local_ipcore,
            &lsms.global.iprint,lsms.global.istop,32);
  }

}
