%NProcShared=4 
%Mem=4GB  
%Chk=H2O_H2O_3080.chk 

#p CCSD(T)/CEP-31G 

H2O_H2O_3080

0 1
O             -0.000000         -0.009833          0.000000
H             -0.799571         -0.580965          0.000000
H              0.799571         -0.580965          0.000000
O             -0.003161          4.867637          0.000000
H              0.570393          5.665473          0.000000
H             -0.943883          5.151445          0.000000

