%NProcShared=4 
%Mem=4GB  
%Chk=O2_O2_163.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

O2_O2_163 

0 3
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      2.234457     2.234457     0.000000 
O      3.054701     3.054701     0.000000 


