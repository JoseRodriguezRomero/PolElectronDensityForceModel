%NProcShared=4 
%Mem=4GB  
%Chk=O2_O2_16.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

O2_O2_16 

0 3
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      0.848528     0.848528     0.000000 
O      1.668772     1.668772     0.000000 


