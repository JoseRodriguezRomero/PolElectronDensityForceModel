#!/bin/bash


FOLDER="OO Coord/"
for METHOD in UFF GAFF Ghemical MMFF94 MMFF94s
do
    echo " " >> $FOLDER$METHOD_energies.txt
    rm -f "$FOLDER""$METHOD""_energies.txt"

    PRINT_VAR=$(obenergy -ff $METHOD "$FOLDER"H2O.xyz | grep "TOTAL ENERGY")
    echo "H2O            "$PRINT_VAR >> $FOLDER$METHOD"_energies.txt"

    for ((i=0; i<5000; i=i+1))
    do
        PRINT_VAR=$(obenergy -ff $METHOD "$FOLDER"H2O_H2O_"$i".xyz | grep "TOTAL ENERGY")
        echo "H2O_H2O_$i     "$PRINT_VAR >> $FOLDER$METHOD"_energies.txt"
    done
    echo "\n" >> $FOLDER$METHOD"_energies.txt"
done

