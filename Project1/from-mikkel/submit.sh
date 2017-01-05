#!/bin/sh

#PBS -q hpcintro
cd $PBS_O_WORKDIR

#echo 'CPU info'
#lscpu

sizes=( 16 23 32 45 64 90 128 181 256 362 512 724 1024 1448 2048 2896 4096 )

for size in ${sizes[@]}
 do
  if [ "$perm" = "blk" ]; then
   ./matmult_c.studio $perm $size $size $size $blksize
  else
   ./matmult_c.studio $perm $size $size $size
  fi
 done
