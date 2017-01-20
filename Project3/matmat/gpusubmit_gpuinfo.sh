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
#BSUB -W 00:01
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

/appl/cuda/8.0/samples/1_Utilities/deviceQuery/deviceQuery

/appl/cuda/8.0/samples/1_Utilities/bandwidthTest/bandwidthTest
