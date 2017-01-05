module load studio
make clean
make OPT='-g -fast' CXX=sunCC

qsub -v perm='mkn' submit.sh -o mkn_o.txt -e mkn_e.txt
qsub -v perm='nkm' submit.sh -o nkm_o.txt -e nkm_e.txt
qsub -v perm='nmk' submit.sh -o nmk_o.txt -e nmk_e.txt
qsub -v perm='mnk' submit.sh -o mnk_o.txt -e mnk_e.txt
qsub -v perm='kmn' submit.sh -o kmn_o.txt -e kmn_e.txt
qsub -v perm='knm' submit.sh -o knm_o.txt -e knm_e.txt

qsub -v perm='lib' submit.sh -o lib_o.txt -e lib_e.txt

qsub -v perm='blk',blksize=40 submit.sh -o blk_o.txt -e blk_e.txt

