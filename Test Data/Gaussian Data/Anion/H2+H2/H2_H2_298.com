%NProcShared=4 
%Mem=4GB  
%Chk=H2_H2_298.chk 

#p ROCCSD(T)/cc-pVTZ 

H2_H2_298 

-1 2
H     -0.740000     0.000000     0.000000 
H      0.000000     0.000000     0.000000 
H      3.507250     3.507250     0.000000 
H      4.030509     4.030509     0.000000 


