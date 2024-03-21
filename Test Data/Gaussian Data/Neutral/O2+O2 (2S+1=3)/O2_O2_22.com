%NProcShared=4 
%Mem=4GB  
%Chk=O2_O2_22.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

O2_O2_22 

0 3
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      0.905097     0.905097     0.000000 
O      1.725341     1.725341     0.000000 


