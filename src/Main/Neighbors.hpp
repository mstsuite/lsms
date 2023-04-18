/* Neighbors - cell-based neighbor listing classes and iterators.
 * Written by David M. Rogers
 * Copyright (C) UT-Battelle, 2022-2023
 *
 * This source code is part of FastParticleToolkit.
 * It has been released under the same BSD-3 open-source license as the rest of
 * LSMS.
 *
 */
#ifndef _NEIGHBORS_HPP
#define _NEIGHBORS_HPP
#include <math.h>
#include <stdint.h>

#include <iterator>

#include "buildLIZandCommLists.hpp"

/** The three integers indexing a lattice point.
 *  {i,j,k} are a union with n[3].
 */
union LatticePt {
  struct { int i,j,k; };
  int n[3];
  LatticePt(int _i, int _j, int _k) : n{_i,_j,_k} {}
  bool operator== (const LatticePt& b) const {
    return i == b.i && j == b.j && k == b.k;
  };
  bool operator!= (const LatticePt& b) const {
    return i != b.i || j != b.j || k != b.k;
  };
};


/** The geometry associated with a Grid of nx,ny,nz cells.
 *
 *  The lattice vectors of the grid are:
 *
 *   * x: Lx
 *   * y: Lyx, Ly
 *   * z: Lzx, Lzy, Lz
 *
 * The lattice vectors indexing cells within the Grid are:
 *   * i: Lx/nx
 *   * j: Lyx/ny, Ly/ny
 *   * k: Lzx/nz, Lzy/nz, Lz/nz
 *
 * Each cell is an axis-aligned rectangle whose lengths
 * are the diagonals of the cell lattice vectors, above.
 *
 */
struct Geom {
  static const int XX = 0;
  static const int YY = 1;
  static const int ZZ = 2;
  static const int YX = 3;
  static const int ZX = 4;
  static const int ZY = 5;
  const Real L[6];
  const Real h[6];
  const int n[3];

  Geom(Real Lx, Real Ly, Real Lz, int nx, int ny, int nz, Real Lyx = 0.0,
       Real Lzx = 0.0, Real Lzy = 0.0)
      : L{Lx, Ly, Lz, Lyx, Lzx, Lzy},
        h{Lx / nx, Ly / ny, Lz / nz, Lyx / ny, Lzx / nz, Lzy / nz},
        n{nx, ny, nz} {}

  inline int cells() const { return n[0] * n[1] * n[2]; }
  void offset(const LatticePt &pt, Real x[3]) const {
    x[0] = h[XX] * pt.i + h[YX] * pt.j + h[ZX] * pt.k;
    x[1] = h[YY] * pt.j + h[ZY] * pt.k;
    x[2] = h[ZZ] * pt.k;
  }

  inline void decodeBin(const unsigned int bin, int &i, int &j, int &k) const {
    i = bin % n[0];
    j = (bin / n[0]) % n[1];
    k = bin / (n[0] * n[1]);
  }

  inline uint32_t calcBin(const int i, const int j, const int k) const {
    return (k * n[1] + j) * n[0] + i;
  }

  inline uint32_t calcBinF(Real &x, Real &y, Real &z) const {
    // 1. wrap into range [0,Lz)
    const auto wz = floor(z / L[ZZ]);
    z -= wz * L[ZZ];
    y -= wz * L[ZY];
    x -= wz * L[ZX];
    const int nz = z / h[ZZ];

    // 2. wrap into range y_0(nz) + [0,Ly)
    const Real y0 = nz * h[ZY];
    const auto wy = floor((y - y0) / L[YY]);
    y -= wy * L[YY];
    x -= wy * L[YX];
    const int ny = (y - y0) / h[YY];

    // 3. wrap into range x_0(ny,nz) + [0,Lx)
    const Real x0 = nz * h[ZX] + ny * h[ZY];
    x -= L[XX] * floor((x - x0) / L[XX]);
    const int nx = (x - x0) / h[XX];

    return calcBin(nx, ny, nz);
  }
};

/** NeighborCells iterates over LatticePt-s within
 *  a distance Rc of a cell at the origin.
 *
 *  This is used during pair computations to iterate
 *  over neighboring cells as in:
 *
 *      NeighborCells nbr(geo, 4.0);
 *      for(const LatticePt &n : nbr) {
 *          // use n.i, n.j, n.k;
 *      }
 */
struct NeighborCells {
  static constexpr Real eps = 1e-6;
  const Geom g;
  const Real Rc;
  const int k0, k1;

  NeighborCells(const Geom &geo, const Real _Rc)
      : g(geo),
        Rc(_Rc),
        k0(ceil(-(geo.h[Geom::ZZ] + Rc) / geo.h[Geom::ZZ] + eps)),
        k1(floor((geo.h[Geom::ZZ] + Rc) / geo.h[Geom::ZZ] - eps)) {}

  class iterator {
    LatticePt n;  // the returned lattice point
    int k1, j1, i1;
    Real R2z;

   public:
    const Geom g;
    const Real Rc;

    using iterator_category = std::forward_iterator_tag;
    using difference_type = int;

    using value_type = LatticePt;
    using reference = LatticePt &;
    using pointer = LatticePt *;

    // ctors:
    iterator(const Geom &geo, const Real _Rc, int _k0, int _k1)
        : g(geo), Rc(_Rc), n(0, 0, _k0 - 1), k1(_k1), j1(0), i1(0) {
      next();
    }
    // iterator(const CellRange *prng) : rng(*prng) {}

    reference operator*() { return n; }
    pointer operator->() { return &n; }
    // Prefix increment
    iterator &operator++() {
      next();
      return *this;
    }
    // call this function until it returns true
    void next() {
      if (++n.i <= i1) return;
      while (true) {
        if (++n.j > j1) {
          if (++n.k > k1) {
            // sets to NeighborCells.end()
            n.i = 0;
            n.j = 0;
            return;
          }
          // find new j-range for this k value
          // Real dz = fabs(n.k*g.h[Geom::ZZ]) - g.h[Geom::YY];
          // dz -= dz*(dz < 0.0f); // max(0, dy) // dist. from z-box edge
          int absk = n.k > 0 ? n.k - 1 : (n.k < 0 ? -n.k - 1 : 0);
          Real dz = absk * g.h[Geom::ZZ];  // equiv. to above, fewer operations
          R2z = Rc * Rc - dz * dz;
          if(R2z < 0.0f) continue; // go to next kReal

          Real y0 = -n.k * g.h[Geom::ZY];
          Real ywid = g.h[Geom::YY] + sqrt(R2z);

          n.j = ceil((y0 - ywid) / g.h[Geom::YY] + NeighborCells::eps);
          j1 = floor((y0 + ywid) / g.h[Geom::YY] - NeighborCells::eps);

        }
        // find new i-range for this j value
        Real dy =
            fabs(n.j * g.h[Geom::YY] + n.k * g.h[Geom::ZY]) - g.h[Geom::YY];
        dy -= dy * (dy < 0.0f);  // max(0, dy) = dist. from y-box edge
        Real R2y = R2z - dy * dy;
        if (R2y < 0.0f)
          continue;
          // goto next j
        Real x0 =
            -n.j * g.h[Geom::YX] - n.k * g.h[Geom::ZX];  // shifted base-pt.
        Real xwid = g.h[Geom::XX] + sqrt(R2y);

        n.i = ceil((x0 - xwid) / g.h[Geom::XX] + NeighborCells::eps);
        i1 = floor((x0 + xwid) / g.h[Geom::XX] - NeighborCells::eps);
        return;
      }
    }
    // Postfix increment
    iterator operator++(int) {
      iterator tmp = *this;
      ++(*this);
      return tmp;
    }

    friend bool operator==(const iterator &a, const iterator &b) {
      return a.n == b.n;
    };
    friend bool operator!=(const iterator &a, const iterator &b) {
      return a.n != b.n;
    };
  };

  iterator begin() { return iterator(g, Rc, k0, k1); }
  iterator end() { return iterator(g, Rc, k1 + 1, k1); }
};

std::vector<std::vector<int>> sort_atoms(CrystalParameters &crystal,
                                         const Geom &geo);
int buildLIZ(CrystalParameters &crystal, int idx, std::vector<LIZInfoType> &LIZ,
             const Geom &geo, const std::vector<std::vector<int>> &cell);
#endif
