#!/bin/sh
#PBS -q hpcintro
#PBS -l nodes=1:ppn=1
#PBS -l walltime=60:00
cd $PBS_O_WORKDIR

# collect -h dch,on,dcm,on,l2h,on,l2m,on ./ass2_main_jacobi
# // ./ass2_main <method type> <NN> <d> <kmax>

N=100
kmax=1000000
method=jacobi
d=0.01

sizesn=( 10 30 50 100 150 200 300 )

for N in ${sizesn[@]}
 do
  ./ass2_main $method $N $d $kmax
done

method=gauss
for N in ${sizesn[@]}
 do
  ./ass2_main $method $N $d $kmax
done

