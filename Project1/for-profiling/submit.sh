#!/bin/sh

#PBS -q hpcintro
#PBS -l nodes=1:ppn=1
#PBS -l walltime=60:00
#PBS -N collector
cd $PBS_O_WORKDIR

module load studio

#echo 'CPU info'
#lscpu

sizes=( 256 1024 )

MFLOPS_MAX_IT=1

for size in ${sizes[@]}
 do
  if [ "$perm" = "blk" ]; then
   collect -h dch,on,dcm,on,l2h,on,l2m,on ./matmult_c.studio $perm $size $size $size $blksize
  else
   collect -h dch,on,dcm,on,l2h,on,l2m,on ./matmult_c.studio $perm $size $size $size
  fi
 done
