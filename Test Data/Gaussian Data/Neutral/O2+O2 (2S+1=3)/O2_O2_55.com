%NProcShared=4 
%Mem=4GB  
%Chk=O2_O2_55.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

O2_O2_55 

0 3
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      1.216224     1.216224     0.000000 
O      2.036468     2.036468     0.000000 


