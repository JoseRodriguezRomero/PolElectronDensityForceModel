%NProcShared=4 
%Mem=4GB  
%Chk=N2_N2_253.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

N2_N2_253 

-1 2
N     -1.090000     0.000000     0.000000 
N      0.000000     0.000000     0.000000 
N      3.082986     3.082986     0.000000 
N      3.853732     3.853732     0.000000 


