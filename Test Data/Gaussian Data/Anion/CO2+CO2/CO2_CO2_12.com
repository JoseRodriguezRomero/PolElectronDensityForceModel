%NProcShared=8 
%Mem=6GB  
%Chk=CO2_CO2_12.chk 

#p ROCCSD(T)/cc-pVTZ 

CO2_CO2_12 

-1 2
O     -2.326000     0.000000     0.000000 
C     -1.163000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      1.160000     0.000000     0.000000 
C      2.323000     0.000000     0.000000 
O      3.486000     0.000000     0.000000 


