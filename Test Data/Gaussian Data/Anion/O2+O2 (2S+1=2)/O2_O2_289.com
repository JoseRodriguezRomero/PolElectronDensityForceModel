%NProcShared=8 
%Mem=6GB  
%Chk=O2_O2_289.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

O2_O2_289 

-1 2
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      3.422397     3.422397     0.000000 
O      4.242641     4.242641     0.000000 


