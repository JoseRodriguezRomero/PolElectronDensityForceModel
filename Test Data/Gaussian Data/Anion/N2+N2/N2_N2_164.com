%NProcShared=4 
%Mem=4GB  
%Chk=N2_N2_164.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

N2_N2_164 

-1 2
N     -1.090000     0.000000     0.000000 
N      0.000000     0.000000     0.000000 
N      0.000000     3.160000     0.000000 
N      0.000000     4.250000     0.000000 

