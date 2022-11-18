/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#ifndef LSMS_LSMSALGORITHMS_H
#define LSMS_LSMSALGORITHMS_H

namespace LSMSAlgorithms {

  const size_t numAlgorithms = 2;

  ///////// lizConstruction /////////
  const size_t lizConstruction = 0;
  const signed char lizConstruction_original = 0;
  const signed char lizConstruction_bricks = 1;
  
  const signed char lizConstruction_default = lizConstruction_original;

  ///////// potentialIO /////////
  const size_t potentialIO = 1;
  const signed char potentialIO_serial = 0;
  const signed char potentialIO_parallel = 1;

  const signed char potentialIO_default = potentialIO_serial;
}

#endif
