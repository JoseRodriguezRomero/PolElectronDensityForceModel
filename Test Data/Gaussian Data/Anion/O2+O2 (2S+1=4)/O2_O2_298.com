%NProcShared=8 
%Mem=6GB  
%Chk=O2_O2_298.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

O2_O2_298 

-1 4
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      3.507250     3.507250     0.000000 
O      4.327494     4.327494     0.000000 

