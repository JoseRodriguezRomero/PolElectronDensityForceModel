%NProcShared=4 
%Mem=4GB  
%Chk=N2_N2_64.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

N2_N2_64 

-1 2
N     -1.090000     0.000000     0.000000 
N      0.000000     0.000000     0.000000 
N      1.301076     1.301076     0.000000 
N      2.071823     2.071823     0.000000 


