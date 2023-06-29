#include "Complex.hpp"
#include <vector>
#include "Array3d.hpp"
#include "Matrix.hpp"
#include "Indices.hpp"

class ClebschGordan {
public:
  ClebschGordan(){};
  static AngularMomentumIndices a;
  static int lmax_cg, kmax_cg;
  static Matrix <Complex> evec;
  static Array3d <Real> cg;
  static Array3d <Real> l3list;
  static void init(int lmax);
  static int callsize(int l1, int m1, int l2, int m2);
  static int kdelta(int k1, int k2);
  static Real calProductSeries(int start_val, int end_val);
  static Real calfx(Real x, int l1, int m1, int l2, int m2);
  static std::vector <Real> calculateClebschGordanCoefficient(int l_1, int m_1, int l_2, int m_2);
  static void populateClebschGordanTable();
  static Real getClebschGordanCoefficient(int l1, int m1, int l2, int m2, int l3, int m3);
};
