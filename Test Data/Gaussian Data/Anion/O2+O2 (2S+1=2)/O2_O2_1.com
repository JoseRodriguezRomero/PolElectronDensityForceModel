%NProcShared=8 
%Mem=6GB  
%Chk=O2_O2_1.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

O2_O2_1 

-1 2
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      0.707107     0.707107     0.000000 
O      1.527351     1.527351     0.000000 


