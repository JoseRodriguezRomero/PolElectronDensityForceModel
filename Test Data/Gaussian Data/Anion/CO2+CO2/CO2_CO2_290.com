%NProcShared=8 
%Mem=6GB  
%Chk=CO2_CO2_290.chk 

#p ROCCSD(T)/cc-pVTZ 

CO2_CO2_290 

-1 2
O     -2.326000     0.000000     0.000000 
C     -1.163000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      0.000000     4.840000     0.000000 
C      0.000000     6.003000     0.000000 
O      0.000000     7.166000     0.000000 

