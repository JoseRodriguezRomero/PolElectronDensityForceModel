%NProcShared=4 
%Mem=4GB  
%Chk=N2_66.chk 

#p CCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)

N2_66 

0 1
N      0.000000     0.000000     0.000000 
N      2.780000     0.000000     0.000000 


