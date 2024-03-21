%NProcShared=4 
%Mem=4GB  
%Chk=O2_O2_46.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

O2_O2_46 

0 3
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      1.131371     1.131371     0.000000 
O      1.951615     1.951615     0.000000 


