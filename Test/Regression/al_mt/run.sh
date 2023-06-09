#!/bin/bash


cwd=$(pwd)

for i in {1..9}
do
   cd ${cwd}/${i}
   mpirun -np 1 ${LSMS_BASE} i_lsms.lua
done
