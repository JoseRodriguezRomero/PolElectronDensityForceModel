%NProcShared=4 
%Mem=4GB  
%Chk=N2_N2_7.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

N2_N2_7 

-1 2
N     -1.090000     0.000000     0.000000 
N      0.000000     0.000000     0.000000 
N      0.763675     0.763675     0.000000 
N      1.534422     1.534422     0.000000 


