%NProcShared=8 
%Mem=6GB  
%Chk=O2_O2_197.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

O2_O2_197 

-1 2
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      0.000000     3.600000     0.000000 
O      0.000000     4.760000     0.000000 


