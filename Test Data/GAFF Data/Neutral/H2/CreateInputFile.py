import copy
import numpy as np

def WriteXYZFile(file_name,dist):
    # Write the XYZ data to file_name based on the specified interatomic H-H
    # bond distance.
    
    H1X = 0.0;
    H1Y = 0.0;
    H1Z = 0.0;
    
    H2X = dist;
    H2Y = 0.0;
    H2Z = 0.0;
    
    fileID = open(file_name,"w+")
    fileID.write("2\n\n")
    
    fileID.write("H  ")
    fileID.write("{:12.6f} ".format(H1X))
    fileID.write("{:12.6f} ".format(H1Y))
    fileID.write("{:12.6f} ".format(H1Z)+"\n")
    
    fileID.write("H  ")
    fileID.write("{:12.6f} ".format(H2X))
    fileID.write("{:12.6f} ".format(H2Y))
    fileID.write("{:12.6f} ".format(H2Z)+"\n")
    
    fileID.write("\n\n")
    fileID.close()
    return

file_num = 0
base_name = "H2_"

for i in range(101):
    dist = 0.8 + i*3.0/100.0
    file_num = file_num + 1
    WriteXYZFile(base_name+str(file_num-1)+".xyz",dist)
