%NProcShared=4 
%Mem=4GB  
%Chk=N2_N2_44.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

N2_N2_44 

-1 2
N     -1.090000     0.000000     0.000000 
N      0.000000     0.000000     0.000000 
N      0.000000     1.560000     0.000000 
N      0.000000     2.650000     0.000000 


