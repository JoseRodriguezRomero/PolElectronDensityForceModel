%NProcShared=4 
%Mem=4GB  
%Chk=H2_H2_64.chk 

#p ROCCSD(T)/cc-pVTZ 

H2_H2_64 

-1 2
H     -0.740000     0.000000     0.000000 
H      0.000000     0.000000     0.000000 
H      1.301076     1.301076     0.000000 
H      1.824335     1.824335     0.000000 


