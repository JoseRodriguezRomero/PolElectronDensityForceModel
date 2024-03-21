%NProcShared=4 
%Mem=4GB  
%Chk=O2_O2_29.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

O2_O2_29 

0 3
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      0.000000     1.360000     0.000000 
O      0.000000     2.520000     0.000000 


