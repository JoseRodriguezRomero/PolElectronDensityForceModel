#!/bin/bash

cd $SCRATCH/Anion/O2+O2\ \(2S+1\=2\)/

#$ -N "O2+O2(2)"
#$ -q short.q
#$ -cwd
#$ -o $1.log
#$ -e $1.err
#$ -l h_rt=10:00:00
#$ -l h_vmem=1.0G
#$ -pe openmp 8
#$ -M jose.rodrigues-romero@student.uibk.ac.at
#$ -m e

export OMP_NUM_THREADS=$NSLOTS
module load gaussian

for i in {0..302}
do
	g16 O2_O2_$i.com
done

