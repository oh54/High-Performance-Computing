#!/bin/sh
module load studio
make clean_all
make CFLAGS='-fast -xopenmp'

# sequential experiments, varying d and varying n
#qsub submit-seqd.sh -o seqd_o.txt -e seqd_e.txt
#qsub submit-seqn.sh -o seqn_o.txt -e seqn_e.txt

# parallel experiments, varying n and varying t
qsub submit-parn.sh -o parn_o.txt -e parn_e.txt
qsub submit-part.sh -o part_o.txt -e part_e.txt
