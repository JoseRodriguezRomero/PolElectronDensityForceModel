%NProcShared=8 
%Mem=6GB  
%Chk=O2_O2_49.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

O2_O2_49 

-1 2
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      1.159655     1.159655     0.000000 
O      1.979899     1.979899     0.000000 

