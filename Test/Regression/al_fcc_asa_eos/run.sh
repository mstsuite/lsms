#!/bin/bash

export LSMS_BASE=/home/fmoitzi/Development/CLionProjects/lsms/cmake-build-debug-wsl/bin/lsms
cwd=$(pwd)

for i in {0..6}; do
  cd ${i}
  rm -rf k.out
  cp w_al.0 v_al.0
  mpirun -np 1 ${LSMS_BASE} i_lsms.lua
  cd ${cwd}
done
