#!/bin/sh
#BSUB -J Bocanegra-FP-mpi
#BSUB -o mpi_output_file
#BSUB -e mpi_error_file
#BSUB -n 20
#BSUB -q ht-10g
#BSUB cwd /home/bocanegra.c/FP/
work=/home/bocanegra.c/FP/
export OMP_NUM_THREADS=60
export OMP_STACKSIZE=10G
export OMP_NESTED=TRUE

cd $work
tempfile1=hostlistrun
tempfile2=hostlist-tcp
echo $LSB_MCPU_HOSTS > $tempfile1
declare -a hosts
read -a hosts < ${tempfile1}
for ((i=0; i<${#hosts[@]}; i += 2)) ;
  do
   HOST=${hosts[$i]}
   CORE=${hosts[(($i+1))]}
   echo $HOST:$CORE >> $tempfile2
done
#####################################################
#####
###DO NOT EDIT ANYTHING ABOVE THIS LINE#########
#####################################################

mpirun -np 15 -prot -TCP -lsf ./main_mpi