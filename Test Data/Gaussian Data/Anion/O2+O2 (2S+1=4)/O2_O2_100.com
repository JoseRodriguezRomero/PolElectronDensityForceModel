%NProcShared=8 
%Mem=6GB  
%Chk=O2_O2_100.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

O2_O2_100 

-1 4
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      1.640488     1.640488     0.000000 
O      2.460732     2.460732     0.000000 

