%NProcShared=4 
%Mem=4GB  
%Chk=N2_79.chk 

#p CCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

N2_79 

0 1
N      0.000000     0.000000     0.000000 
N      3.170000     0.000000     0.000000 


