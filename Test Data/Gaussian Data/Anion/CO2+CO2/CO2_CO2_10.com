%NProcShared=8 
%Mem=6GB  
%Chk=CO2_CO2_10.chk 

#p ROCCSD(T)/cc-pVTZ 

CO2_CO2_10 

-1 2
O     -2.326000     0.000000     0.000000 
C     -1.163000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      0.791960     0.791960     0.000000 
C      1.614325     1.614325     0.000000 
O      2.436690     2.436690     0.000000 

