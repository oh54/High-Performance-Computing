#!/bin/sh
#PBS -q hpcintro
#PBS -l nodes=1:ppn=1
#PBS -l walltime=60:00
cd $PBS_O_WORKDIR

# collect -h dch,on,dcm,on,l2h,on,l2m,on ./ass2_main_jacobi
# // ./ass2_main <method type> <NN> <d> <kmax>

N=100
kmax=100000
method=jacobi

sizesd=( 0.1 0.01 0.001 0.0001 0.00001 0.000001 0.0000001 0.00000001 0.000000001 0.0000000001 0.00000000001 )

for d in ${sizesd[@]}
 do
  ./ass2_main $method $N $d $kmax
done

method=gauss
for d in ${sizesd[@]}
 do
  ./ass2_main $method $N $d $kmax
done

