#!/bin/bash

for i in {66..302}
do
	psi4 CO2_CO2_$i.inp -n 11
done

