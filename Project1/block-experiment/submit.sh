#!/bin/sh

#PBS -q hpcintro
#PBS -l nodes=1:ppn=1
#PBS -l walltime=60:00

cd $PBS_O_WORKDIR

module load studio

#echo 'CPU info'
#lscpu

sizes=( 32 64 512 2048 )
block_sizes=( 2 20 40 60 100 300 )

for size in ${sizes[@]}
 do
  for blksize in ${block_sizes[@]}
   do
    ./matmult_c.studio blk $size $size $size $blksize
   done
 done
