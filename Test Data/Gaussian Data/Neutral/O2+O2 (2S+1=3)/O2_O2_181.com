%NProcShared=4 
%Mem=4GB  
%Chk=O2_O2_181.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

O2_O2_181 

0 3
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      2.404163     2.404163     0.000000 
O      3.224407     3.224407     0.000000 


