%NProcShared=4 
%Mem=4GB  
%Chk=O2_O2_15.chk 

#p ROCCSD(T)/cc-pVTZ 

O2_O2_15 

0 5
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      1.200000     0.000000     0.000000 
O      2.360000     0.000000     0.000000 


