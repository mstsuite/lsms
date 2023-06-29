#include <iostream>
#include "ClebschGordan.hpp"

int ClebschGordan::lmax_cg;
int ClebschGordan::kmax_cg;
Matrix <Complex> ClebschGordan::evec;
Array3d <Real> ClebschGordan::cg,ClebschGordan::l3list;
AngularMomentumIndices ClebschGordan::a;
void ClebschGordan::init(int lmax){
//  Calculates Clebsch-Gordan coefficient <l1m1l2m2|lm>
//  This is done using a recursive algorithm

    const Complex sqrtm1(0.0, 1.0);
//  Setting the three vectors ep1, em1 and e0
//  ep1 = -(x + iy)/sqrt(2)
//  em1 = +(x - iy)/sqrt(2)
//  e0 = z
//  Matrix evec : evec(0,:) = em1, evec(1:,) = e0, evec(2:,) = ep1
    evec.resize(3, 3);
    evec(0,0) = 1.0/sqrt(2);
    evec(0,1) = -sqrtm1/sqrt(2);
    evec(0,2) = 0.0;
    evec(1,0) = 0.0;
    evec(1,1) = 0.0;
    evec(1,2) = 1.0;  
    evec(2,0) = -1.0/sqrt(2);
    evec(2,1) = -sqrtm1/sqrt(2);
    evec(2,2) = 0.0;

//  Allocating memory for Clebsch-Gordan storage clb
//  and |l1 - l2| <= l3 <= l1 + l2 + 1 
    lmax_cg = lmax;
    kmax_cg = (lmax+1)*(lmax+1);
    cg.resize(kmax_cg, kmax_cg, 2*lmax_cg + 2);
    l3list.resize(kmax_cg, kmax_cg, 2*lmax_cg + 2);

    ClebschGordan::a.init(lmax_cg);
    ClebschGordan::populateClebschGordanTable();
}

int ClebschGordan::callsize(int l1, int m1, int l2, int m2){
   int lsize = 0;
   if (abs(m1 + m2) <= abs(l1 - l2)) {
     lsize = 2*fmin(l1,l2) + 2;
   } 
   else {
     lsize = l1 + l2 - abs(m1 + m2) + 2;
   }
   return lsize;
}

int ClebschGordan::kdelta(int k1, int k2){
   if (k1 == k2){
     return 1;
   }
   else {
     return 0;
   }
}

Real ClebschGordan::calProductSeries(int start_val, int end_val){
//  Returns start_val!/end_val!
    Real prod = 1.0;
    for(int j=start_val;j>end_val;j=j-1){
       prod = prod*j;
    }
    return prod;
}

Real ClebschGordan::calfx(Real x, int l1, int m1, int l2, int m2){
    Real pre_A, A;
    int m = m1 + m2;
    pre_A = (pow(x,2) - pow(m,2))*(pow(l1 + l2 + 1, 2) - pow(x, 2))*
          (pow(x, 2) - pow(l1 - l2, 2))/(pow(x, 2) * (4*pow(x, 2) - 1.0));
    A = sqrt(pre_A);
    return A;
}

std::vector <Real> ClebschGordan::calculateClebschGordanCoefficient(int l1, int m1, int l2, int m2){
    int temp;
    Real A_ml, A_pl, A_ol;
    std::vector <Real> cgl;
    std::vector <int> l3l;
    int lsize = ClebschGordan::callsize(l1, m1, l2, m2);
    cgl.resize(lsize);
    l3l.resize(lsize);

//  <l1l2|l1+l2+1> = 0
    cgl[lsize-1] = 0.0;
    
    l3l[lsize-1] = l1 + l2 + 1;
    for(int i=lsize-2;i>=0;i=i-1){
       l3l[i] = l3l[i+1] - 1;
    }

    Real coeff = ClebschGordan::calProductSeries(2*l1, l1)/
           ClebschGordan::calProductSeries(2*l1 + 2*l2, 2*l2);
    coeff = coeff*ClebschGordan::calProductSeries(l1, 0);
    
    if ((l1 - m1) % 2 == 0){
      coeff = coeff*(ClebschGordan::calProductSeries(l1 + l2 + m1 + m2, l1 + m1)/
                     ClebschGordan::calProductSeries(l1 - m1, (l1 - m1)/2));
      coeff = coeff*(1.0/ClebschGordan::calProductSeries((l1 - m1)/2, 0));    
    }
    else {
      coeff = coeff*(ClebschGordan::calProductSeries(l1 + l2 + m1 + m2, l1 + m1)/
                     ClebschGordan::calProductSeries(l1 - m1, (l1 - m1 + 1)/2));
      coeff = coeff*(1.0/ClebschGordan::calProductSeries((l1 - m1 + 1)/2, 0));
    }

    if ((l2 + m2) % 2 == 0){
      coeff = coeff*(ClebschGordan::calProductSeries(l1 + l2 - m1 - m2, l2 - m2)/
                     ClebschGordan::calProductSeries(l2 + m2, (l2 + m2)/2));
      coeff = coeff*(1.0/ClebschGordan::calProductSeries((l2 + m2)/2, 0));
    }
    else {
      coeff = coeff*(ClebschGordan::calProductSeries(l1 + l2 - m1 - m2, l2 - m2)/
                     ClebschGordan::calProductSeries(l2 + m2, (l2 + m2 + 1)/2));
      coeff = coeff*(1.0/ClebschGordan::calProductSeries((l2 + m2 + 1)/2, 0));
    }


//  coeff = <l1m1l2m2|l1+l2,m1+m2>
    cgl[lsize-2] = sqrt(coeff);

    for(int i=lsize-2;i>=1;i=i-1){
       temp = l3l[i];
       A_ol = m1 - m2 + (m1 + m2)*(l2*(l2 + 1.0) - l1*(l1 + 1.0))/(temp*(temp + 1.0));
       A_pl = ClebschGordan::calfx(temp + 1, l1, m1, l2, m2);
       A_ml = ClebschGordan::calfx(temp, l1, m1, l2, m2);
       cgl[i - 1] = (A_ol/A_ml)*cgl[i] - (A_pl/A_ml)*cgl[i + 1];
    }
    return cgl;
}

void ClebschGordan::populateClebschGordanTable(){
    int k1,k2,l1,l2,m1,m2,lsize,i;
    std::vector <Real> cgtemp;    

    for(k1=0;k1<kmax_cg;k1++){
       for(k2=0;k2<kmax_cg;k2++){
          l1 = ClebschGordan::a.lofk[k1];
          l2 = ClebschGordan::a.lofk[k2];

          m1 = ClebschGordan::a.mofk[k1];
          m2 = ClebschGordan::a.mofk[k2];
          lsize = ClebschGordan::callsize(l1, m1, l2, m2);
          cgtemp = ClebschGordan::calculateClebschGordanCoefficient(l1, m1, l2, m2);
          for(i=0;i<lsize;i++){
             cg(k1, k2, i) = cgtemp[i];
          }
          l3list(k1, k2, lsize-1) = l1 + l2 + 1;
          for(i=lsize-2;i>=0;i=i-1){
             l3list(k1, k2, i) = l3list(k1, k2, i+1) - 1;
          }
       }
    }
}

Real ClebschGordan::getClebschGordanCoefficient(int l1, int m1, int l2, int m2, int l3, int m3){
    int k1, k2, lsize,i;
    Real cgcoeff = 0.0;
    if (m1 + m2 != m3){
      return 0;
    }
    else if (abs(m1) > l1 || abs(m2) > l2 || abs(m3) > l3){
      return 0;
    }
    else {
      k1 = pow(l1, 2) + l1 + m1;
      k2 = pow(l2, 2) + l2 + m2;
      lsize = ClebschGordan::callsize(l1, m1, l2, m2);
      for (i=0;i<lsize;i++){
         if (l3 == l3list(k1, k2, i)){
           cgcoeff = cg(k1, k2, i);
         }
      }
      return cgcoeff;
    }
}
