%NProcShared=4 
%Mem=4GB  
%Chk=H2_H2_99.chk 

#p ROCCSD(T)/cc-pVTZ 

H2_H2_99 

-1 2
H     -0.740000     0.000000     0.000000 
H      0.000000     0.000000     0.000000 
H      2.320000     0.000000     0.000000 
H      3.060000     0.000000     0.000000 


