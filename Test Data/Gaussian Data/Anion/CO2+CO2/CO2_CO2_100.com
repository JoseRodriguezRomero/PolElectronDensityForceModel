%NProcShared=8 
%Mem=6GB  
%Chk=CO2_CO2_100.chk 

#p ROCCSD(T)/cc-pVTZ 

CO2_CO2_100 

-1 2
O     -2.326000     0.000000     0.000000 
C     -1.163000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      1.640488     1.640488     0.000000 
C      2.462853     2.462853     0.000000 
O      3.285218     3.285218     0.000000 

