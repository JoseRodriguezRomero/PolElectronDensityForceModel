%NProcShared=4 
%Mem=4GB  
%Chk=O2_O2_19.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

O2_O2_19 

0 3
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      0.876812     0.876812     0.000000 
O      1.697056     1.697056     0.000000 


