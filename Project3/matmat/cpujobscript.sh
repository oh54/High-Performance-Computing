#!/bin/sh
### General options
### –- specify queue --
#BSUB -q hpc
### -- set the job Name --
#BSUB -J My_Application
### -- ask for number of cores (default: 1) --
#BSUB -n 4
### -- set walltime limit: hh:mm --
#BSUB -W 00:10
### -- set the email address --
# please uncomment the following line and put in your e-mail address,
# if you want to receive e-mail notifications on a non-default address
##BSUB -u your_email_address
### -- send notification at start --
#BSUB -B
### -- send notification at completion --
#BSUB -N
### -- Specify the output and error file. %J is the job-id --
### -- -o and -e mean append, -oo and -eo mean overwrite --
#BSUB -o Output_%J.out
#BSUB -e Error_%J.err
#BSUB -R "select[model=XeonE5_2660v3]"

# here follow the commands you want to execute
module load cuda/8.0

sizes=( 16 32 48 64 96 128 176 256 352 512 704 1024 1408 2048 2816 4096 )

for size in ${sizes[@]}
 do
    ./matmult_f.nvcc2 lib $size $size $size
 done


