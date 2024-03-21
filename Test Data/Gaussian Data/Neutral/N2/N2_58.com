%NProcShared=4 
%Mem=4GB  
%Chk=N2_58.chk 

#p CCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

N2_58 

0 1
N      0.000000     0.000000     0.000000 
N      2.540000     0.000000     0.000000 


