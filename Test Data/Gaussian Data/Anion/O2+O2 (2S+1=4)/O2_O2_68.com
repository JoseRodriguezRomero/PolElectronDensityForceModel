%NProcShared=8 
%Mem=6GB  
%Chk=O2_O2_68.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

O2_O2_68 

-1 4
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      0.000000     1.880000     0.000000 
O      0.000000     3.040000     0.000000 

