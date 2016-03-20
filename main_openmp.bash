#!/bin/sh
#BSUB -J Bocanegra‐hello
#BSUB -o output_file
#BSUB -e error_file
#BSUB -n 20
#BSUB -q ht-10g
#BSUB cwd /home/bocanegra.c/FP/
work=/home/bocanegra.c/FP/
cd $work
export OMP_NUM_THREADS=20
export OMP_STACKSIZE=10G
./main_openmp
