#!/bin/sh
#PBS -q hpcintro
#PBS -l nodes=1:ppn=4
#PBS -l walltime=60:00
cd $PBS_O_WORKDIR

# collect -h dch,on,dcm,on,l2h,on,l2m,on ./ass2_main_jacobi
# // ./ass2_main <method type> <NN> <d> <kmax>
# threads

kmax=1000000
method=omp
d=0.01

sizesn=( 10 30 50 100 150 200 300 500 1000 2000 5000 10000 )

for N in ${sizesn[@]}
 do
  OMP_NUM_THREADS=4 ./ass2_main $method $N $d $kmax
done

method=jacobi
for N in ${sizesn[@]}
 do
  OMP_NUM_THREADS=1 ./ass2_main $method $N $d $kmax
done

method=omp2
for N in ${sizesn[@]}
 do
  ./ass2_main $method $N $d $kmax
done

