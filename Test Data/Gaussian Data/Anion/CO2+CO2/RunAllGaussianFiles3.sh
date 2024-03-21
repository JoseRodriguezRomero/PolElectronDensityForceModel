#!/bin/bash

cd $SCRATCH/Anion/CO2+CO2

#$ -N "CO2+CO2"
#$ -q short.q
#$ -cwd
#$ -o $1.log
#$ -e $1.err
#$ -l h_rt=10:00:00
#$ -l h_vmem=1.5G
#$ -pe openmp 8
#$ -M jose.rodrigues-romero@student.uibk.ac.at
#$ -m e

export OMP_NUM_THREADS=$NSLOTS
module load gaussian

for i in {301..302}
do
	g16 CO2_CO2_$i.com
done



