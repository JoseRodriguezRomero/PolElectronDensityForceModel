%NProcShared=4 
%Mem=4GB  
%Chk=O2_O2_76.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

O2_O2_76 

0 3
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      1.414214     1.414214     0.000000 
O      2.234457     2.234457     0.000000 


