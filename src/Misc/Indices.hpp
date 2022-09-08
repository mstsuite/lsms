#ifndef LSMS_INDICES_HPP
#define LSMS_INDICES_HPP

#include <vector>

namespace lsms {

template <typename Integer,
          typename = std::enable_if_t<std::is_integral<Integer>::value>>
Integer get_jmax(Integer lmax) {
  return (lmax + 1) * (lmax + 2) / 2;
}

template <typename Integer,
          typename = std::enable_if_t<std::is_integral<Integer>::value>>
Integer get_kmax(Integer lmax) {
  return (lmax + 1) * (lmax + 1);
}

}  // namespace lsms

class AngularMomentumIndices {
 public:
  static int lmax, ndlj, ndlm;
  static std::vector<int> lofk, mofk, lofj, mofj, kofj;

  static void init(int _lmax) {
    lmax = _lmax;
    ndlj = (lmax + 1) * (lmax + 1);
    ndlm = ((lmax + 1) * (lmax + 2)) / 2;

    mofk.resize(ndlj);
    lofk.resize(ndlj);

    mofj.resize(ndlm);
    lofj.resize(ndlm);
    kofj.resize(ndlm);

    int j = 0;
    int k = 0;

    for (int l = 0; l <= lmax; l++) {
      for (int m = 0; m <= l; m++) {
        lofj[j] = l;
        mofj[j] = m;
        kofj[j] = (l + 1) * (l + 1) - l + m - 1;

        j++;
      }

      for (int m = -l; m <= l; m++) {
        lofk[k] = l;
        mofk[k] = m;

        // printf("k,l,m: %d %d %d, lofk,mofk:%d %d\n",k,l,m,lofk[k],mofk[k]);

        k++;
      }
    }
    // printf("k=%d\n",k);
  }
};

#endif
