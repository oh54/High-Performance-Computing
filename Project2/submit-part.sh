#!/bin/sh
#PBS -q hpcintro
#PBS -l nodes=1:ppn=20
#PBS -l walltime=60:00
cd $PBS_O_WORKDIR

# collect -h dch,on,dcm,on,l2h,on,l2m,on ./ass2_main_jacobi
# // ./ass2_main <method type> <NN> <d> <kmax>
# threads

kmax=1000000
method=omp
d=0.01
N=10000
sizest=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 )

for t in ${sizest[@]}
 do
  OMP_NUM_THREADS=$t ./ass2_main $method $N $d $kmax
done

#method=omp2
#for t in ${sizest[@]}
# do
#  OMP_NUM_THREADS=$t ./ass2_main $method $N $d $kmax
#done
