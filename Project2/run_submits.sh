#!/bin/sh
module load studio
make clean_all
make

qsub submit.sh -o jacobi_o.txt -e jacobi_e.txt
