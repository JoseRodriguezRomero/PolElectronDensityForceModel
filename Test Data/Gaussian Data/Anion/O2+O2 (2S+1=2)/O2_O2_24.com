%NProcShared=8 
%Mem=6GB  
%Chk=O2_O2_24.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

O2_O2_24 

-1 2
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      1.320000     0.000000     0.000000 
O      2.480000     0.000000     0.000000 


