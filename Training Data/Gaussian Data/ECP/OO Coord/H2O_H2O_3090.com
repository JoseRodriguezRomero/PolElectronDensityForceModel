%NProcShared=4 
%Mem=4GB  
%Chk=H2O_H2O_3090.chk 

#p CCSD(T)/CEP-31G 

H2O_H2O_3090

0 1
O             -0.000000         -0.009833          0.000000
H             -0.799571         -0.580965          0.000000
H              0.799571         -0.580965          0.000000
O             -0.008176          4.863788          0.000000
H             -0.038837          5.845912          0.000000
H             -0.927273          4.516273          0.000000
