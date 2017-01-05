#!/bin/sh

#PBS -q hpcintro
cd $PBS_O_WORKDIR


echo 'CPU info'
lscpu

#echo 'Running studio NKM 50 60 70'
#./matmult_c.studio nkm 50 60 70

#sizes=( 16 32 64 128 256 512 1024 2048 )
sizes=( 16 23 32 45 64 90 128 181 256 362 512 724 1024 1448 2048 )

for size in ${sizes[@]}
 do
  ./matmult_c.studio $perm $size $size $size
 done

#for N in ${sizes[@]}
# do 
#  for K in ${sizes[@]}
#   do
#    for M in ${sizes[@]}
#     do
#      echo $M $N $K
#      ./matmult_c.studio mnk $M $N $K
#     done
#   done
# done
