%NProcShared=4 
%Mem=4GB  
%Chk=O2_O2_12.chk 

#p ROCCSD(T)/cc-pVTZ 

O2_O2_12 

0 5
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      1.160000     0.000000     0.000000 
O      2.320000     0.000000     0.000000 


