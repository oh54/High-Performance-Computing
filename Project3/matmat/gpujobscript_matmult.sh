#!/bin/sh
### General options
### –- specify queue --
#BSUB -q k40
### -- set the job Name --
#BSUB -J k40job
### -- ask for number of cores (default: 1) --
#BSUB -n 2
### -- Select the resources: 2 gpus in exclusive process mode --
#BSUB -R "rusage[ngpus_excl_p=1]"
### -- set walltime limit: hh:mm --
#BSUB -W 00:10
### -- set the email address --
# please uncomment the following line and put in your e-mail address,
# if you want to receive e-mail notifications on a non-default address
##BSUB -u your_email_address
### -- send notification at start --
#BSUB -B
### -- send notification at completion--
#BSUB -N
### -- Specify the output and error file. %J is the job-id --
### -- -o and -e mean append, -oo and -eo mean overwrite --
#BSUB -o gpu-%J.out
#BSUB -e gpu_%J.err

# -- end of LSF options --

# Load the cuda module
module load cuda/8.0

# * this program sees only the gpu's which are requested.
# * double check the number of your requested gpus above.
# * max-walltime in this queue is 30 minutes
# * this node has only 12 cores, so please don't request more
#   than 2 cpu-cores
#
sizes1=( 16 32 48 64 96 128 176 256 )

for size in ${sizes1[@]}
 do
    ./matmult_f.nvcc2 gpu1 $size $size $size
 done

sizes=( 16 32 48 64 96 128 176 256 352 512 704 1024 1408 2048 2816 4096 )

for size in ${sizes[@]}
 do
    ./matmult_f.nvcc2 gpu2 $size $size $size
 done

for size in ${sizes[@]}
 do
    ./matmult_f.nvcc2 gpu3 $size $size $size
 done

for size in ${sizes[@]}
 do
    ./matmult_f.nvcc2 gpu4 $size $size $size
 done

for size in ${sizes[@]}
 do
    ./matmult_f.nvcc2 gpu5 $size $size $size
 done

for size in ${sizes[@]}
 do
    ./matmult_f.nvcc2 gpulib $size $size $size
 done
#./myprogram
