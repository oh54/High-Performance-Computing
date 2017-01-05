module load studio
make clean
make OPT='-g -fast' CXX=sunCC

qsub -v perm='mkn' submit.sh -o mkn_o.txt -e mkn_e.txt
qsub -v perm='nkm' submit.sh -o nkm_o.txt -e nkm_e.txt
