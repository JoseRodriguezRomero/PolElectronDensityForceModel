%NProcShared=4 
%Mem=4GB  
%Chk=O2_O2_61.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

O2_O2_61 

0 3
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      1.272792     1.272792     0.000000 
O      2.093036     2.093036     0.000000 


