%NProcShared=4 
%Mem=4GB  
%Chk=H2_H2_4.chk 

#p ROCCSD(T)/cc-pVTZ 

H2_H2_4 

-1 2
H     -0.740000     0.000000     0.000000 
H      0.000000     0.000000     0.000000 
H      0.735391     0.735391     0.000000 
H      1.258650     1.258650     0.000000 


