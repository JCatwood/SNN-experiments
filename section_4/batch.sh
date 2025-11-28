#!/bin/bash

for scene in 1 2 3; do
    Rscript SNN_m_cmp_unknown.R 1 $scene 0 1>>tmp.out 2>>tmp.out
done
