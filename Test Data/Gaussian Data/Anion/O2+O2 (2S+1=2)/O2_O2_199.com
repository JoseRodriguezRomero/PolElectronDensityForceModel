%NProcShared=8 
%Mem=6GB  
%Chk=O2_O2_199.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

O2_O2_199 

-1 2
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      2.573869     2.573869     0.000000 
O      3.394113     3.394113     0.000000 


