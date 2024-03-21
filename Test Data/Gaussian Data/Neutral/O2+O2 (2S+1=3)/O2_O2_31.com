%NProcShared=4 
%Mem=4GB  
%Chk=O2_O2_31.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

O2_O2_31 

0 3
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      0.989949     0.989949     0.000000 
O      1.810193     1.810193     0.000000 


