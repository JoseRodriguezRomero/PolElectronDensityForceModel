%NProcShared=8 
%Mem=6GB  
%Chk=O2_O2_274.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

O2_O2_274 

-1 2
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      3.280975     3.280975     0.000000 
O      4.101219     4.101219     0.000000 


