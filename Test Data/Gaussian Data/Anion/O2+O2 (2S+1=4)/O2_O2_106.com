%NProcShared=8 
%Mem=6GB  
%Chk=O2_O2_106.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

O2_O2_106 

-1 4
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      1.697056     1.697056     0.000000 
O      2.517300     2.517300     0.000000 


