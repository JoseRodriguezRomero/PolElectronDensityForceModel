%NProcShared=4 
%Mem=4GB  
%Chk=O2_O2_298.chk 

#p ROCCSD(T)/cc-pVTZ 

O2_O2_298 

0 5
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      3.507250     3.507250     0.000000 
O      4.327494     4.327494     0.000000 


