%NProcShared=4 
%Mem=4GB  
%Chk=H2_H2_1.chk 

#p CCSD(T)/cc-pVTZ 

H2_H2_1 

-1 2
H     -0.740000     0.000000     0.000000 
H      0.000000     0.000000     0.000000 
H      0.707107     0.707107     0.000000 
H      1.230366     1.230366     0.000000 

