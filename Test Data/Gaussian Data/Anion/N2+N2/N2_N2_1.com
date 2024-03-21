%NProcShared=4 
%Mem=4GB  
%Chk=N2_N2_1.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

N2_N2_1 

-1 2
N     -1.090000     0.000000     0.000000 
N      0.000000     0.000000     0.000000 
N      0.707107     0.707107     0.000000 
N      1.477853     1.477853     0.000000 


