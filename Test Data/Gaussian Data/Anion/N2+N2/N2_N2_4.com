%NProcShared=4 
%Mem=4GB  
%Chk=N2_N2_4.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

N2_N2_4 

-1 2
N     -1.090000     0.000000     0.000000 
N      0.000000     0.000000     0.000000 
N      0.735391     0.735391     0.000000 
N      1.506137     1.506137     0.000000 

