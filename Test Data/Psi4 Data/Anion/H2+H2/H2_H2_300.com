%NProcShared=4 
%Mem=4GB  
%Chk=H2_H2_300.chk 

#p CCSD(T)/cc-pVTZ 

H2_H2_300 

-1 2
H     -0.740000     0.000000     0.000000 
H      0.000000     0.000000     0.000000 
H      5.000000     0.000000     0.000000 
H      5.740000     0.000000     0.000000 


