#!/usr/bin/env bash


for i in {0..3}
do
   echo "w_al.$i" "v_al.$i"
   cp "w_al.$i" "v_al.$i"
done