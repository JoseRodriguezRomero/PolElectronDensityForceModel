%NProcShared=4 
%Mem=4GB  
%Chk=N2_N2_229.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

N2_N2_229 

-1 2
N     -1.090000     0.000000     0.000000 
N      0.000000     0.000000     0.000000 
N      2.856711     2.856711     0.000000 
N      3.627458     3.627458     0.000000 


