%NProcShared=4 
%Mem=4GB  
%Chk=O2_O2_265.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

O2_O2_265 

0 3
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      3.196123     3.196123     0.000000 
O      4.016367     4.016367     0.000000 


