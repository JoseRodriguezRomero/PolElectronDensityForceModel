%NProcShared=4 
%Mem=4GB  
%Chk=N2_N2_205.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

N2_N2_205 

-1 2
N     -1.090000     0.000000     0.000000 
N      0.000000     0.000000     0.000000 
N      2.630437     2.630437     0.000000 
N      3.401184     3.401184     0.000000 

