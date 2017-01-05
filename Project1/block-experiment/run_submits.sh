module load studio
make clean
make OPT='-g -fast' CXX=sunCC

qsub submit.sh -o blk_o.txt -e blk_e.txt

