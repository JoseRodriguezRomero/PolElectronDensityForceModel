%NProcShared=8 
%Mem=6GB  
%Chk=O2_O2_190.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

O2_O2_190 

-1 4
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      2.489016     2.489016     0.000000 
O      3.309260     3.309260     0.000000 


