%NProcShared=4 
%Mem=4GB  
%Chk=N2_N2_214.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

N2_N2_214 

-1 2
N     -1.090000     0.000000     0.000000 
N      0.000000     0.000000     0.000000 
N      2.715290     2.715290     0.000000 
N      3.486036     3.486036     0.000000 

