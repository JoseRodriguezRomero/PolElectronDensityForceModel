%NProcShared=8 
%Mem=6GB  
%Chk=O2_O2_157.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

O2_O2_157 

-1 2
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      2.177889     2.177889     0.000000 
O      2.998133     2.998133     0.000000 


