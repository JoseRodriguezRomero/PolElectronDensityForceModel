%NProcShared=8 
%Mem=6GB  
%Chk=O2_O2_7.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

O2_O2_7 

-1 4
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      0.763675     0.763675     0.000000 
O      1.583919     1.583919     0.000000 


