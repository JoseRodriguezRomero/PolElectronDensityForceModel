%NProcShared=4 
%Mem=4GB  
%Chk=O2_O2_256.chk 

#p ROCCSD(T)/cc-pVTZ 

O2_O2_256 

0 5
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      3.111270     3.111270     0.000000 
O      3.931514     3.931514     0.000000 

