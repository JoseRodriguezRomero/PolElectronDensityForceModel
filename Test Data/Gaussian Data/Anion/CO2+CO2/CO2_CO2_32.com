%NProcShared=8 
%Mem=6GB  
%Chk=CO2_CO2_32.chk 

#p ROCCSD(T)/cc-pVTZ 

CO2_CO2_32 

-1 2
O     -2.326000     0.000000     0.000000 
C     -1.163000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      0.000000     1.400000     0.000000 
C      0.000000     2.563000     0.000000 
O      0.000000     3.726000     0.000000 


