%NProcShared=4 
%Mem=4GB  
%Chk=O2_O2_28.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

O2_O2_28 

0 3
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      0.961665     0.961665     0.000000 
O      1.781909     1.781909     0.000000 


