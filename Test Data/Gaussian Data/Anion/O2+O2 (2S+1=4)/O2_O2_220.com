%NProcShared=8 
%Mem=6GB  
%Chk=O2_O2_220.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

O2_O2_220 

-1 4
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      2.771859     2.771859     0.000000 
O      3.592102     3.592102     0.000000 


