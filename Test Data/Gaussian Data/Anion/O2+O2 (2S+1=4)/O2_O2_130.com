%NProcShared=8 
%Mem=6GB  
%Chk=O2_O2_130.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

O2_O2_130 

-1 4
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      1.923330     1.923330     0.000000 
O      2.743574     2.743574     0.000000 

