#!/bin/bash

touch EFP_energies.txt
rm EFP_energies.txt
touch EFP_energies.txt

for ((i=0; i<5000; i++))
do
    PRINT_VAR=$(efpmd "OH Coord/H2O_H2O_$i.inp"  | grep "TOTAL ENERGY")
    echo "$PRINT_VAR" >> EFP_energies.txt
done


for ((i=0; i<5000; i++))
do
    PRINT_VAR=$(efpmd "OO Coord/H2O_H2O_$i.inp" | grep "TOTAL ENERGY")
    echo "$PRINT_VAR" >> EFP_energies.txt
done

