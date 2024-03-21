#!/bin/bash

MB_POL_DIR=~/Downloads/MBX-master/bin/

touch MB_Pol_energies.txt
rm MB_Pol_energies.txt
touch MB_Pol_energies.txt

PRINT_VAR=$("$MB_POL_DIR"single_point "H2O.nrg")
echo "$PRINT_VAR" >> MB_Pol_energies.txt

for ((i=0; i<5000; i++))
do
    PRINT_VAR=$("$MB_POL_DIR"single_point "OH Coord/H2O_H2O_$i.nrg")
    echo "$PRINT_VAR" >> MB_Pol_energies.txt
done


for ((i=0; i<5000; i++))
do
    PRINT_VAR=$("$MB_POL_DIR"single_point "OO Coord/H2O_H2O_$i.nrg")
    echo "$PRINT_VAR" >> MB_Pol_energies.txt
done

