%NProcShared=4 
%Mem=4GB  
%Chk=N2_N2_157.chk 

#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

N2_N2_157 

-1 2
N     -1.090000     0.000000     0.000000 
N      0.000000     0.000000     0.000000 
N      2.177889     2.177889     0.000000 
N      2.948635     2.948635     0.000000 

