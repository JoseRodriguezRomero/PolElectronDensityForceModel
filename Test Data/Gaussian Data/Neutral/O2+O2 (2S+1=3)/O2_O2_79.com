%NProcShared=4 
%Mem=4GB  
%Chk=O2_O2_79.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

O2_O2_79 

0 3
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      1.442498     1.442498     0.000000 
O      2.262742     2.262742     0.000000 


