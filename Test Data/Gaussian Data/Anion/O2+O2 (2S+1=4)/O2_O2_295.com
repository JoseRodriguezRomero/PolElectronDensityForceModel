%NProcShared=8 
%Mem=6GB  
%Chk=O2_O2_295.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

O2_O2_295 

-1 4
O     -1.160000     0.000000     0.000000 
O      0.000000     0.000000     0.000000 
O      3.478965     3.478965     0.000000 
O      4.299209     4.299209     0.000000 


