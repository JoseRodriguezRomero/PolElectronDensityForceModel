%NProcShared=4 
%Mem=4GB  
%Chk=O2_O2_196.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

O2_O2_196 

0 3
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      2.545584     2.545584     0.000000 
O      3.365828     3.365828     0.000000 


