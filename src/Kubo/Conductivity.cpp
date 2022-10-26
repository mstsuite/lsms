#include "Conductivity.hpp"

Conductivity::Conductivity(LSMSSystemParameters &lsms, LSMSCommunication &comm, LocalTypeInfo &local, Real volume) {
    omega = volume;
    pola = lsms.n_spin_pola;
    efermi = Complex(lsms.chempot,0.0001);
    cm.resize(local.atom.size(), pola);
    sigmatilde1.resize(3, 3, pola);
    sigmatilde2.resize(3, 3, pola);
    sigmatilde3.resize(3, 3, pola);
    sigmatilde4.resize(3, 3, pola);
    sigma.resize(3, 3, pola);
    rho.resize(3, 3, pola);
    spin_summed_rho.resize(3, 3);

    double timeSingleScatterers=MPI_Wtime();
    local.tmatStore=0.0;
    local.JxStore=0.0;
    local.JyStore=0.0;
    local.JzStore=0.0;
    expectTmatCommunication(comm,local);
    expectJxCommunication(comm,local);
    expectJyCommunication(comm,local);
    expectJzCommunication(comm,local); 
    
    if(lsms.global.iprint>=1) printf("calculate current and t matrix.\n");
     
    for (int i=0;i<local.atom.size();i++){
       for(int is=0;is<pola;is++){
          cm(i,is).init(lsms, local, local.atom[i], i, efermi, is);
       }
    }

    if(lsms.global.iprint>=2) printf("About to send t matrices\n");
     sendTmats(comm,local);
    if(lsms.global.iprint>=2) printf("About to send J matrices\n");
     sendJx(comm,local);
     sendJy(comm,local);
     sendJz(comm,local);
    if(lsms.global.iprint>=2) printf("About to finalize t matrices communication\n");
      finalizeTmatCommunication(comm);
    if(lsms.global.iprint>=1) printf("Received all t matrices");
    if(lsms.global.iprint>=2) printf("About to finalize J matrices communication\n");
      finalizeJxCommunication(comm);
      finalizeJyCommunication(comm);
      finalizeJzCommunication(comm);
    if(lsms.global.iprint>=2) printf("Recieved all J matricies\n");
    timeSingleScatterers=MPI_Wtime()-timeSingleScatterers;
    if(lsms.global.iprint>=0) printf("timeSingleScatteres = %lf sec\n",timeSingleScatterers);
    /*
    if (comm.rank == 0){
    std::cout << std::endl;
    std::cout << "Local Atom on Rank" << comm.rank << " : Current Mat" << std::endl;
    std::cout << std::endl;
    for (int i=0;i<16;i++){
       for (int j=0;j<16;j++){
          std::cout << local.tmatStore(j+i*16,0) << "   ";
       }
       std::cout << std::endl;
    }   
    std::cout << std::endl;
    std::cout << "Remote Atom on Rank" << comm.rank << " : Current Mat" << std::endl;
    std::cout << std::endl;
    for (int i=0;i<16;i++){
       for (int j=0;j<16;j++){
          std::cout << local.tmatStore(j+i*16,1) << "   ";
       }
       std::cout << std::endl;
    }   
    std::cout << std::endl;
    }
    */
   double timeTauMatrix=MPI_Wtime();
   for (int i=0;i<local.atom.size();i++){
      for (int is=0;is<pola;is++){
        cm(i,is).calTauFull(lsms,local,local.atom[i]);
      }
   }
   
   timeTauMatrix=MPI_Wtime()-timeTauMatrix;
   if(lsms.global.iprint>=0) printf("timeTauMatrix = %lf sec\n",timeTauMatrix);
   fflush(stdout);

   double timeConductivity=MPI_Wtime();
   calSigma(comm,local);
   timeConductivity=MPI_Wtime()-timeConductivity;
   if (comm.rank == 0) {
     writeSigmaTildeMat(sigmatilde1, "sigmatilde1");
     writeSigmaTildeMat(sigmatilde2, "sigmatilde2");
     writeSigmaTildeMat(sigmatilde3, "sigmatilde3");
     writeSigmaTildeMat(sigmatilde4, "sigmatilde4");
     writeRhoMat();
   }
   if(lsms.global.iprint>=0) printf("timeConductivity = %lf sec\n",timeConductivity);
}

void Conductivity::processCurrentMatrix(LocalTypeInfo &local, int index, int dir, Matrix <Complex> &Jout1, int kkrsz, int etype){
    if (etype == 1) {
      if (dir == 1){
        for(int i=0;i<kkrsz;i++){
          for(int j=0;j<kkrsz;j++){
            Jout1(i,j) = local.JxStore(i+j*kkrsz,index);
          }
        }
      }
      else if (dir == 2){
        for(int i=0;i<kkrsz;i++){
          for(int j=0;j<kkrsz;j++){
            Jout1(i,j) = local.JyStore(i+j*kkrsz,index);
          }
        }
      }
      else if (dir == 3){
        for(int i=0;i<kkrsz;i++){
          for(int j=0;j<kkrsz;j++){
            Jout1(i,j) = local.JzStore(i+j*kkrsz,index);
          }
        }
      }
    }
    else if (etype == 2) {
      if (dir == 1){
        for(int i=0;i<kkrsz;i++){
          for(int j=0;j<kkrsz;j++){
            Jout1(i,j) = pow(-1.0, ClebschGordan::a.lofk[i])*local.JxStore(i+j*kkrsz,index);
          }
        }
      }
      else if (dir == 2){
        for(int i=0;i<kkrsz;i++){
          for(int j=0;j<kkrsz;j++){
            Jout1(i,j) = pow(-1.0, ClebschGordan::a.lofk[i])*local.JyStore(i+j*kkrsz,index);
          }
        }
      }
      else if (dir == 3){
        for(int i=0;i<kkrsz;i++){
          for(int j=0;j<kkrsz;j++){
            Jout1(i,j) = pow(-1.0, ClebschGordan::a.lofk[i])*local.JzStore(i+j*kkrsz,index);
          }
        }
      }
    }
    else if (etype == 3) {
      if (dir == 1) {
        for(int i=0;i<kkrsz;i++){
          for(int j=0;j<kkrsz;j++){
            Jout1(i,j) = pow(-1.0, ClebschGordan::a.lofk[j])*local.JxStore(i+j*kkrsz,index);
          }
        }
      }
      else if (dir == 2){
        for(int i=0;i<kkrsz;i++){
          for(int j=0;j<kkrsz;j++){
            Jout1(i,j) = pow(-1.0, ClebschGordan::a.lofk[j])*local.JyStore(i+j*kkrsz,index);
          }
        }
      }
      else if (dir == 3){
        for(int i=0;i<kkrsz;i++){
          for(int j=0;j<kkrsz;j++){
            Jout1(i,j) = pow(-1.0, ClebschGordan::a.lofk[j])*local.JzStore(i+j*kkrsz,index);
          }
        }
      }
    }
    else if (etype == 4) {
      if (dir == 1){
        for(int i=0;i<kkrsz;i++){
          for(int j=0;j<kkrsz;j++){
            Jout1(i,j) = pow(-1.0, ClebschGordan::a.lofk[i]+ClebschGordan::a.lofk[j])*local.JxStore(i+j*kkrsz,index);
          }
        }
      }
      else if (dir == 2) {
        for(int i=0;i<kkrsz;i++){
          for(int j=0;j<kkrsz;j++){
            Jout1(i,j) = pow(-1.0, ClebschGordan::a.lofk[i]+ClebschGordan::a.lofk[j])*local.JyStore(i+j*kkrsz,index);
          }
        }
      }
      else if (dir == 3) {
        for(int i=0;i<kkrsz;i++){
          for(int j=0;j<kkrsz;j++){
            Jout1(i,j) = pow(-1.0, ClebschGordan::a.lofk[i]+ClebschGordan::a.lofk[j])*local.JzStore(i+j*kkrsz,index);
          }
        }
      }
    }
}

void Conductivity::processTauMatrix(Matrix <Complex> &tau1, Matrix <Complex> &tau2, int m, int n, int is, int kkrsz, int etype){
    if (etype == 1) {
      for (int i=0;i<kkrsz;i++){
        for (int j=0;j<kkrsz;j++){
           tau1(i,j) = cm(m,is).tau1(i,n*kkrsz+j);
           tau2(i,j) = cm(m,is).tau1(n*kkrsz+i,j);
        }
      }
    }
    else if (etype == 2) {
      for (int i=0;i<kkrsz;i++){
        for (int j=0;j<kkrsz;j++){
          tau1(i,j) = cm(m,is).tau1(i,n*kkrsz+j);
          tau2(i,j) = pow(-1.0, ClebschGordan::a.lofk[i] - ClebschGordan::a.lofk[j])*
                      conj(cm(m,is).tau1(j,n*kkrsz+i));
        }
      }
    }
    else if (etype == 3) {
       for (int i=0;i<kkrsz;i++){
         for (int j=0;j<kkrsz;j++){
            tau1(i,j) = pow(-1.0, ClebschGordan::a.lofk[i] - ClebschGordan::a.lofk[j])*
                        conj(cm(m,is).tau1(n*kkrsz+j,i));
            tau2(i,j) = cm(m,is).tau1(n*kkrsz+i,j);
         }
       }
    }
    else if (etype == 4) {
       for (int i=0;i<kkrsz;i++){
         for (int j=0;j<kkrsz;j++){
            tau1(i,j) = pow(-1.0, ClebschGordan::a.lofk[i] - ClebschGordan::a.lofk[j])*
                        conj(cm(m,is).tau1(n*kkrsz+j,i));
            tau2(i,j) = pow(-1.0, ClebschGordan::a.lofk[i] - ClebschGordan::a.lofk[j])*
                        conj(cm(m,is).tau1(j,n*kkrsz+i));
         }
       }
    }
}

Complex Conductivity::calSigmaTilde(LocalTypeInfo &local, int dir1, int dir2, int is, int etype){
    Matrix <Complex> J1(local.atom[0].kkrsz, local.atom[0].kkrsz);
    Matrix <Complex> J2(local.atom[0].kkrsz, local.atom[0].kkrsz);
    Matrix <Complex> taumn(local.atom[0].kkrsz, local.atom[0].kkrsz);
    Matrix <Complex> taunm(local.atom[0].kkrsz, local.atom[0].kkrsz);

    Complex sigmatilde = 0.0;
    for(int m=0;m<local.atom.size();m++){
       if (etype == 2 || etype == 3) {
         processCurrentMatrix(local,m,dir1,J1,local.atom[m].kkrsz,5-etype);
       }
       else {
         processCurrentMatrix(local,m,dir1,J1,local.atom[m].kkrsz,etype);
       }
       for (int n=0;n<local.atom[m].numLIZ;n++){
          processCurrentMatrix(local,local.atom[m].LIZStoreIdx[n],dir2,J2,local.atom[m].kkrsz,etype);
          processTauMatrix(taumn, taunm, m, n, is,local.atom[m].kkrsz, etype);
          for (int L1=0;L1<local.atom[m].kkrsz;L1++){
            for (int L2=0;L2<local.atom[m].kkrsz;L2++){
              for (int L3=0;L3<local.atom[m].kkrsz;L3++){
                for (int L4=0;L4<local.atom[m].kkrsz;L4++){
                  sigmatilde = sigmatilde + J1(L4,L1)*taumn(L1,L2)*J2(L2,L3)*taunm(L3,L4);
                }
              }
            }
          }
       }
    }
    return -sigmatilde/(M_PI*omega);
}

void Conductivity::calSigma(LSMSCommunication &comm, LocalTypeInfo &local){   
   int dirs=3;
   Complex temp; Matrix <Real> spin_summed_sigma(dirs, dirs); 
   Matrix <Real> work(dirs,dirs);
   int ipiv[dirs], info;
   Real sigmatmpsum;

  #pragma omp parallel for collapse(3) private(temp)
   for(int is=0;is<pola;is++){
     for(int dir1=1;dir1<1+dirs;dir1++){
       for(int dir2=1;dir2<1+dirs;dir2++){
         sigmatilde1(dir1-1,dir2-1,is) = calSigmaTilde(local,dir1,dir2,is,1);
         sigmatilde2(dir1-1,dir2-1,is) = calSigmaTilde(local,dir1,dir2,is,2);
         sigmatilde3(dir1-1,dir2-1,is) = calSigmaTilde(local,dir1,dir2,is,3);
         sigmatilde4(dir1-1,dir2-1,is) = calSigmaTilde(local,dir1,dir2,is,4);
         temp = 0.25*(sigmatilde1(dir1-1,dir2-1,is) - sigmatilde2(dir1-1,dir2-1,is)
                  - sigmatilde3(dir1-1,dir2-1,is) + sigmatilde4(dir1-1,dir2-1,is));
         sigma(dir1-1,dir2-1,is) = 0.0230384174*temp.real();
       }
     }
   }

   globalSum(comm, &sigma(0,0,0), 9*pola);

   for(int i=0;i<dirs;i++){
     for(int j=0;j<dirs;j++){
       sigmatmpsum = 0.0;
       for(int is=0;is<pola;is++){
          sigmatmpsum = sigmatmpsum + sigma(i,j,is);
       }
       if (pola == 1) {
          sigmatmpsum = 2*sigmatmpsum;
       }
       spin_summed_sigma(i,j)=sigmatmpsum;
     }
   }
 // Inverting spin-resolved conductivities
   for (int is=0;is<pola;is++){
     invertConductivityMatrix(is);
   }

 // Inverting spin-summed conductivity
   LAPACK::dgetrf_(&dirs, &dirs, &spin_summed_sigma(0, 0), &dirs, &ipiv[0], &info);
   LAPACK::dgetri_(&dirs, &spin_summed_sigma(0,0), &dirs, &ipiv[0], &work(0,0), &dirs, &info);
   for(int i=0;i<dirs;i++){
     for(int j=0;j<dirs;j++){
       spin_summed_rho(i,j) = spin_summed_sigma(i,j);
     }
   }
}

void Conductivity::invertConductivityMatrix(int is){
   int dirs = 3;
   Matrix <Real> tmpinvert(dirs, dirs);
   Matrix <Real> work(dirs, dirs);
   int ipiv[dirs], info;

   for(int i=0;i<dirs;i++){
     for(int j=0;j<dirs;j++){
       tmpinvert(i,j) = sigma(i,j,is);
     }
   }
   LAPACK::dgetrf_(&dirs, &dirs, &tmpinvert(0, 0), &dirs, &ipiv[0], &info);
   LAPACK::dgetri_(&dirs, &tmpinvert(0,0), &dirs, &ipiv[0], &work(0,0), &dirs, &info);
   for (int i=0;i<dirs;i++){
     for (int j=0;j<dirs;j++){
       rho(i,j,is) = tmpinvert(i,j);
     }
   }
}

void Conductivity::writeRhoMat(){
  int dirs=3;
  for (int is=0;is<pola;is++){
    std::cout << "RESISTIVITY SPIN " << is+1 << " (in microOhm-cm)" << std::endl;
    for(int i=0;i<dirs;i++){
      for (int j=0;j<dirs;j++){
        std::cout << rho(i,j,is) << "  ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
  std::cout << "TOTAL RESISTIVITY (in microOhm-cm)" << std::endl;
  for (int i=0;i<dirs;i++){
    for (int j=0;j<dirs;j++){
      std::cout << spin_summed_rho(i,j) << "   ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

void Conductivity::writeSigmaTildeMat(Array3d <Complex> &stm, std::string matname){
  int dirs=3;
  for (int is=0;is<pola;is++){
    std::cout << "Matname Spin " << is+1 << " : " << matname << std::endl;
    for(int i=0;i<dirs;i++){
      for(int j=0;j<dirs;j++){
        std::cout << stm(i,j,is) << "   ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
}
