%NProcShared=4 
%Mem=4GB  
%Chk=O2_O2_67.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

O2_O2_67 

0 3
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      1.329361     1.329361     0.000000 
O      2.149605     2.149605     0.000000 

