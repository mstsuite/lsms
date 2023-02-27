#include "CurrentMatrix.hpp"


void CurrentMatrix::init(LSMSSystemParameters &lsms, LocalTypeInfo &local, 
                         AtomData &a,int lindex, Complex en, int ispin){
    is = ispin;
    local_index = lindex;
    energy = en;
    prel = sqrt(en);
    pnrel = sqrt(en*(1.0 + en*c2inv));
    atom = &a;
    kkrsz = a.kkrsz;
    nrmat = a.nrmat;
    lmax_cg = 2*a.lmax;
    solutionNonRel.init(lsms, a, &local.tmatStore(0,local_index));
    ClebschGordan::init(lmax_cg);

    Jx.retarget(kkrsz, kkrsz, &local.JxStore(0,local_index));
    Jy.retarget(kkrsz, kkrsz, &local.JyStore(0,local_index));
    Jz.retarget(kkrsz, kkrsz, &local.JzStore(0,local_index));
 
    zlrd.resize(a.r_mesh.size(), a.lmax+1, 2);
    calculateSingleScattererSolution(lsms,a,a.vr,energy,prel,pnrel,solutionNonRel);
    calRadialSolutionDerivative(a);
    assembleJxFromRadialIntegral(a);
    calJyzFromJx();

    m.resize(nrmat, nrmat);
    tau0.resize(kkrsz, kkrsz);
    tau1.resize(nrmat, nrmat);
}

Complex CurrentMatrix::calPrefactor(int L, int Lp, int dir, int choice){
    int l, m, lp, mp; 
    Complex coeffsum = 0.0;
    Real cgcoeff = 0.0;
 
    l = ClebschGordan::a.lofk[L];
    m = ClebschGordan::a.mofk[L];

    lp = ClebschGordan::a.lofk[Lp];
    mp = ClebschGordan::a.mofk[Lp];
    if (choice == 0 && lp < lmax_cg) {
      for(int k=-1;k<=1;k++){
        cgcoeff = ClebschGordan::getClebschGordanCoefficient(lp + 1, mp - k, 1, k, lp, mp);
        coeffsum = coeffsum + cgcoeff*ClebschGordan::kdelta(l,lp+1)*
                                      ClebschGordan::kdelta(m,mp-k)*ClebschGordan::evec(k+1,dir-1);
      }
      coeffsum = coeffsum*sqrt((lp + 1.0)/(2*lp + 1.0));
    }
    else if (choice == 1) {
       for(int k=-1;k<=1;k++){
         if (lp == mp && (k == 0 || k == -1)){
           cgcoeff = 0.0;
         } 
         else if (lp == -mp && (k == 0 || k == 1)){
           cgcoeff = 0.0;
         }
         else if (lp == mp + 1 && k == -1) {
           cgcoeff = 0.0;
         }
         else if (lp == -mp + 1 && k == 1) {
           cgcoeff = 0.0;
         }
         else {
           cgcoeff = ClebschGordan::getClebschGordanCoefficient(lp - 1, mp - k, 1, k, lp, mp);
         }
         coeffsum = coeffsum + cgcoeff*ClebschGordan::kdelta(l,lp-1)*
                                       ClebschGordan::kdelta(m,mp-k)*ClebschGordan::evec(k+1,dir-1);
       }
       coeffsum = coeffsum*sqrt((lp*1.0)/(2*lp + 1.0));
    }
    return coeffsum;
}

void CurrentMatrix::calRadialSolutionDerivative(AtomData &a){
    for(int i=0;i<=a.lmax;i++){
          zlrd(0,i,is) = (solutionNonRel.zlr(0,i,is) 
                        - solutionNonRel.zlr(1,i,is))/(a.r_mesh[0] - a.r_mesh[1]);
          for(int j=1;j<a.r_mesh.size()-1;j++){
             zlrd(j,i,is) = (solutionNonRel.zlr(j+1,i,is) 
                           - solutionNonRel.zlr(j-1,i,is))/
                           (a.r_mesh[j+1] - a.r_mesh[j-1]);
          }
          zlrd(a.r_mesh.size()-1,i,is) = (solutionNonRel.zlr(a.r_mesh.size()-1,i,is) 
               - solutionNonRel.zlr(a.r_mesh.size()-2,i,is))/
                (a.r_mesh[a.r_mesh.size()-1] - a.r_mesh[a.r_mesh.size()-2]);
    }
}

Complex CurrentMatrix::calRadialIntegral(AtomData &a, int L, int Lp, int dir, int choice){
    std::vector <Complex> integrand;
    std::vector <Real> integrand_real, integral_real, integrand_imag, integral_imag;
    integrand.resize(a.r_mesh.size());
    integrand_real.resize(a.r_mesh.size());
    integral_real.resize(a.r_mesh.size());
    integrand_imag.resize(a.r_mesh.size());
    integral_imag.resize(a.r_mesh.size());
    Complex integrated(0.0,0.0);
    Complex coeff = calPrefactor(L, Lp, dir, choice);
    int l = ClebschGordan::a.lofk[L];
    int lp = ClebschGordan::a.lofk[Lp];

    if (choice == 0){
      Complex coeff = calPrefactor(L, Lp, dir, 0);
      for(int j=0;j<a.r_mesh.size();j++){
         integrand[j] = coeff*pow(a.r_mesh[j], 2)*solutionNonRel.zlr(j, l, is)*(zlrd(j, lp, is) 
                   - ((1.0*lp)/a.r_mesh[j])*solutionNonRel.zlr(j, lp, is));
      }
    }
    else if (choice == 1) {
      Complex coeff = calPrefactor(L, Lp, dir, 1);
      for(int j=0;j<a.r_mesh.size();j++){
         integrand[j] = coeff*pow(a.r_mesh[j], 2)*solutionNonRel.zlr(j, l, is)*(zlrd(j, lp, is) 
                   + ((lp + 1.0)/a.r_mesh[j])*solutionNonRel.zlr(j, lp, is));
      }
    }
    else {
      Complex coeff(0.0,0.0);
    }

    for(int j=0;j<a.r_mesh.size();j++){
      integrand_real[j] = integrand[j].real();
      integrand_imag[j] = integrand[j].imag();
    }

    integrateOneDim(a.r_mesh, integrand_real, integral_real);
    integrateOneDim(a.r_mesh, integrand_imag, integral_imag);
    integrated = Complex(integral_real[a.jmt-1],integral_imag[a.jmt-1]);
    integrated = integrated -
      0.5*pow(a.r_mesh[a.jmt-1],2)*coeff*solutionNonRel.zlr(a.jmt-1,l,is)*solutionNonRel.zlr(a.jmt-1,lp,is);
    return integrated;
}

void CurrentMatrix::assembleJxFromRadialIntegral(AtomData &a){
    int L1,L2;
    Complex pref(0.0, -2.0*sqrt(2.0));
    for(int L1=0;L1<kkrsz;L1++){
       for(int L2=0;L2<kkrsz;L2++){
          Jx(L1,L2) = pref*(-calRadialIntegral(a,L1,L2,1,0)
                         + calRadialIntegral(a,L1,L2,1,1));
       }
    }
}

void CurrentMatrix::calJyzFromJx(){
    int L,Lp,l,m,lp,mp,lpp,mpp;
    Complex pcoeff;
    for(int L=0;L<kkrsz;L++){
      for(int Lp=0;Lp<kkrsz;Lp++){
         m = ClebschGordan::a.mofk[L];
         mp = ClebschGordan::a.mofk[Lp];
         pcoeff = Complex(0.0,-1.0*(m-mp));
         Jy(L,Lp) = pcoeff*Jx(L,Lp);
      }
    }

    for(int L=0;L<kkrsz;L++){
      for(int Lp=0;Lp<kkrsz;Lp++){
        Jz(L,Lp) = 0.0;
        for(int Lpp=0;Lpp<kkrsz;Lpp++){
          l = ClebschGordan::a.lofk[L];
          m = ClebschGordan::a.mofk[L];
          lp = ClebschGordan::a.lofk[Lp];
          mp = ClebschGordan::a.mofk[Lp];
          lpp = ClebschGordan::a.lofk[Lpp];
          mpp = ClebschGordan::a.mofk[Lpp];
          Jz(L,Lp) = Jz(L,Lp) + Complex(0.0,-0.5)*(sqrt(lpp*(lpp+1) - mpp*(mpp+1))*
            Complex(1.0*ClebschGordan::kdelta(m,mpp+1),0.0) + sqrt(lpp*(lpp+1) - mpp*(mpp-1))*
            Complex(1.0*ClebschGordan::kdelta(m,mpp-1),0.0))*Complex(1.0*ClebschGordan::kdelta(l,lpp),0.0)*Jy(Lpp,Lp)
            - Complex(0.0,-0.5)*(sqrt(lp*(lp+1) - mp*(mp+1))*Complex(1.0*ClebschGordan::kdelta(mpp,mp+1),0.0)
            + sqrt(lp*(lp+1) - mp*(mp-1))*Complex(1.0*ClebschGordan::kdelta(mpp,mp-1),0.0))*
              Complex(1.0*ClebschGordan::kdelta(lpp,lp),0.0)*Jy(L,Lpp);
        }
      }
    }
}

void CurrentMatrix::calTauFull(LSMSSystemParameters &lsms, LocalTypeInfo &local,
                         AtomData &a){

     lsms.angularMomentumIndices.init(2*lsms.maxlmax);
     unsigned int buildKKRMatrixKernel = lsms.global.linearSolver;
     unsigned int linearSolver = lsms.global.linearSolver;
     switch(linearSolver) {
        case MST_LINEAR_SOLVER_ZGETRF:
            buildKKRMatrix(lsms,local,a,is,energy,prel,0,m);
            solveTauFullzgetrf(lsms,local,a,m,tau1,0); break;
        case MST_LINEAR_SOLVER_ZGETRF_CUSOLVER:
 	    deviceStorage->allocate(kkrsz,lsms.n_spin_cant,a.numLIZ,lsms.global.GPUThreads);
	    deviceStorage->allocateAdditional(kkrsz,lsms.n_spin_cant,a.numLIZ,lsms.global.GPUThreads);
	    devM = deviceStorage->getDevM();
	    devT = deviceStorage->getDevTFull();
	    transferFullTMatrixToGPUCUDA(devT, lsms, local, a, is);
	    printf("entering buildKKRMatrixCuda:\n");
            buildKKRMatrixCuda(lsms, local, a, *deviceStorage, deviceAtoms[local_index], is, 0, energy, prel,
                         devM);
	    printf("entering solveTauFullzgetrf_cusolver:\n");
            //solveTauFullzgetrf_cublas(lsms, local, *deviceStorage, a, devT, devM, tau1);
            solveTauFullzgetrf_cusolver(lsms, local, *deviceStorage, a, devT, devM, tau1, is);
        case MST_LINEAR_SOLVER_ZGETRF_CUBLAS:
            deviceStorage->allocate(kkrsz,lsms.n_spin_cant,a.numLIZ,lsms.global.GPUThreads);
            deviceStorage->allocateAdditional(kkrsz,lsms.n_spin_cant,a.numLIZ,lsms.global.GPUThreads);
            devM = deviceStorage->getDevM();
            devT = deviceStorage->getDevTFull();
            transferFullTMatrixToGPUCUDA(devT, lsms, local, a, is);
            printf("entering buildKKRMatrixCuda:\n");
            buildKKRMatrixCuda(lsms, local, a, *deviceStorage, deviceAtoms[local_index], is, 0, energy, prel,
                         devM);
            printf("entering solveTauFullzgetrf_cublas:\n");
            solveTauFullzgetrf_cublas(lsms, local, *deviceStorage, a, devT, devM, tau1);
            //solveTauFullzgetrf_cusolver(lsms, local, *deviceStorage, a, devT, devM, tau1, is);
//        std::cout << "GPU implementation in progress!" << std::endl; break;
     }
     /*
     calculateTauMatrix(lsms, local, a, local_index,
                        is, energy, prel,&tau0(0,0),m,0);
    
    for(int i=0;i<kkrsz;i++){
      for (int j=0;j<kkrsz;j++){
         std::cout << tau0(i,j) << "  ";
      }
      std::cout << std::endl;
    }
    */
    std::cout << "00 block of full matrix" << std::endl;
    std::cout << std::endl;
    for (int i=0; i<kkrsz;i++){
      for (int j=0; j<kkrsz;j++){
         std::cout << tau1(i,j) << "  ";
      }
      std::cout << std::endl;
    }
}

