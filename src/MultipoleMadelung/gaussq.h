//
// Created by F.Moitzi on 15.12.2021.
//

#ifndef MADELUNG_GAUSSQ_H
#define MADELUNG_GAUSSQ_H

extern "C" {
void gaussq(int *kind, int *n, int *kpts, double *x1, double *x2, double *t,
            double *w);
}

#endif  // MADELUNG_GAUSSQ_H
