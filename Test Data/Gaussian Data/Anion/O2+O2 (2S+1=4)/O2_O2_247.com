%NProcShared=8 
%Mem=6GB  
%Chk=O2_O2_247.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

O2_O2_247 

-1 4
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      3.026417     3.026417     0.000000 
O      3.846661     3.846661     0.000000 

