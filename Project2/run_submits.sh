#!/bin/sh
module load studio
make clean_all
make

#qsub submit-seqd.sh -o seqd_o.txt -e seqd_e.txt
qsub submit-seqn.sh -o seqn_o.txt -e seqn_e.txt
