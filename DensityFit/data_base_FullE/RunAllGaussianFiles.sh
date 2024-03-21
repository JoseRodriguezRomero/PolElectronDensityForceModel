#!/bin/bash

module load gaussian
input_files=$(ls | grep .com)

for input_file in $input_files
do
    aux_name=$(echo $input_file | grep -o '^[^\.]*')
    
    g16 $input_file

    echo $aux_name

    formchk "$aux_name".chk "$aux_name".fchk
    cubegen 1 fdensity "$aux_name".fchk "$aux_name"_X.cub -1 h < cubegenX.dat
    cubegen 1 fdensity "$aux_name".fchk "$aux_name"_Y.cub -1 h < cubegenY.dat
    cubegen 1 fdensity "$aux_name".fchk "$aux_name"_Z.cub -1 h < cubegenZ.dat
done

