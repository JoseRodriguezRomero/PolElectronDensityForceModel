%NProcShared=4 
%Mem=4GB  
%Chk=N2_N2_136.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

N2_N2_136 

-1 2
N     -1.090000     0.000000     0.000000 
N      0.000000     0.000000     0.000000 
N      1.979899     1.979899     0.000000 
N      2.750645     2.750645     0.000000 


