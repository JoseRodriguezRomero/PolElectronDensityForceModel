%NProcShared=8 
%Mem=6GB  
%Chk=O2_O2_214.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

O2_O2_214 

-1 4
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      2.715290     2.715290     0.000000 
O      3.535534     3.535534     0.000000 


