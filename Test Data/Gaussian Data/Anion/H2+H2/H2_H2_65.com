%NProcShared=4 
%Mem=4GB  
%Chk=H2_H2_65.chk 

#p ROCCSD(T)/cc-pVTZ 

H2_H2_65 

-1 2
H     -0.740000     0.000000     0.000000 
H      0.000000     0.000000     0.000000 
H      0.000000     1.840000     0.000000 
H      0.000000     2.580000     0.000000 


