%NProcShared=4 
%Mem=4GB  
%Chk=H2_H2_125.chk 

#p ROCCSD(T)/cc-pVTZ 

H2_H2_125 

-1 2
H     -0.740000     0.000000     0.000000 
H      0.000000     0.000000     0.000000 
H      0.000000     2.640000     0.000000 
H      0.000000     3.380000     0.000000 


