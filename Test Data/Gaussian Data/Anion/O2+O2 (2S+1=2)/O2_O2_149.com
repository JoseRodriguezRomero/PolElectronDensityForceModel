%NProcShared=8 
%Mem=6GB  
%Chk=O2_O2_149.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

O2_O2_149 

-1 2
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      0.000000     2.960000     0.000000 
O      0.000000     4.120000     0.000000 


