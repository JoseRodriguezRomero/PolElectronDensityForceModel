%NProcShared=4 
%Mem=4GB  
%Chk=O2_O2_34.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

O2_O2_34 

0 3
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      1.018234     1.018234     0.000000 
O      1.838478     1.838478     0.000000 


