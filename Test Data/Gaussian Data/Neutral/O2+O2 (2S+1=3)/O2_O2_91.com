%NProcShared=4 
%Mem=4GB  
%Chk=O2_O2_91.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

O2_O2_91 

0 3
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      1.555635     1.555635     0.000000 
O      2.375879     2.375879     0.000000 


