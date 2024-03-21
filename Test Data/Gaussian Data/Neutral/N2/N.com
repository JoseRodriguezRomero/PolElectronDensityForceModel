%NProcShared=4 
%Mem=4GB  
%Chk=N2_0.chk 

#p CCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

N 

0 4
N      0.000000     0.000000     0.000000 


