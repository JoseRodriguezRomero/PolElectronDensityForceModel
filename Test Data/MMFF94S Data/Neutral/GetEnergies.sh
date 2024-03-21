#!/bin/bash

for MOLEC in CO2
do
    # MOLEC + MOLEC Energies
    FOLDER="$MOLEC+$MOLEC/"
    METHOD="MMFF94S"
    
    cd $FOLDER
    python CreateInputFile.py
    cd ..
    
    echo " " >> $FOLDER$METHOD_energies.txt
    rm -f "$FOLDER""$METHOD""_energies.txt"
    
    PRINT_VAR=$(obenergy -ff $METHOD "$FOLDER$MOLEC".xyz | grep "TOTAL ENERGY")
    echo "$MOLEC          "$PRINT_VAR >> $FOLDER$METHOD"_energies.txt"
    
    for ((i=0; i<303; i=i+1))
    do
        PRINT_VAR=$(obenergy -ff $METHOD "$FOLDER""$MOLEC"_"$MOLEC"_"$i".xyz | grep "TOTAL ENERGY")
        echo "$MOLEC"_"$MOLEC"_"$i     "$PRINT_VAR >> $FOLDER$METHOD"_energies.txt"
    done
    echo "\n" >> $FOLDER$METHOD"_energies.txt"
    
    ## MOLEC Energies
    FOLDER="$MOLEC/"
    METHOD="MMFF94S"
    
    cd $FOLDER
    python CreateInputFile.py
    cd ..
    
    echo " " >> $FOLDER$METHOD_energies.txt
    rm -f "$FOLDER""$METHOD""_energies.txt"
    
    for ((i=0; i<101; i=i+1))
    do
        PRINT_VAR=$(obenergy -ff $METHOD "$FOLDER""$MOLEC"_"$i".xyz | grep "TOTAL ENERGY")
        echo "$MOLEC"_"$i     "$PRINT_VAR >> $FOLDER$METHOD"_energies.txt"
    done
    echo "\n" >> $FOLDER$METHOD"_energies.txt"
done
