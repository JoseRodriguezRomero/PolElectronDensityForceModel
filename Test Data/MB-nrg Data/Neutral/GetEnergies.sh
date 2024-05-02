#!/bin/bash

for MOLEC in CO2
do
    # MOLEC + MOLEC Energies
    FOLDER="$MOLEC+$MOLEC/"
    METHOD="MB-nrg"
    
    cd $FOLDER
    python CreateInputFile.py
    cd ..
    
    echo " " >> $FOLDER$METHOD_energies.txt
    rm -f "$FOLDER""$METHOD""_energies.txt"
    
    PRINT_VAR=$(/Users/joseantoniorodriguesromero/Downloads/MBX/src/single_point "$FOLDER$MOLEC".nrg | grep "Energy=")
    echo "$MOLEC          "$PRINT_VAR >> $FOLDER$METHOD"_energies.txt"
    
    for ((i=0; i<303; i=i+1))
    do
        PRINT_VAR=$(/Users/joseantoniorodriguesromero/Downloads/MBX/src/single_point "$FOLDER""$MOLEC"_"$MOLEC"_"$i".nrg | grep "Energy=")
        echo "$MOLEC"_"$MOLEC"_"$i     "$PRINT_VAR >> $FOLDER$METHOD"_energies.txt"
    done
    echo "\n" >> $FOLDER$METHOD"_energies.txt"
done
