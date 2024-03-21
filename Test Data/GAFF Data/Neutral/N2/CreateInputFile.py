import copy
import numpy as np

def WriteXYZFile(file_name,dist):
    # Write the XYZ data to file_name based on the specified O...O interatomic
    # O-O bond distance.
    
    N1X = 0.0;
    N1Y = 0.0;
    N1Z = 0.0;
    
    N2X = dist;
    N2Y = 0.0;
    N2Z = 0.0;
    
    fileID = open(file_name,"w+")
    fileID.write("2\n\n")
    
    fileID.write("N  ")
    fileID.write("{:12.6f} ".format(N1X))
    fileID.write("{:12.6f} ".format(N1Y))
    fileID.write("{:12.6f} ".format(N1Z)+"\n")
    
    fileID.write("N  ")
    fileID.write("{:12.6f} ".format(N2X))
    fileID.write("{:12.6f} ".format(N2Y))
    fileID.write("{:12.6f} ".format(N2Z)+"\n")
    
    fileID.write("\n\n")
    fileID.close()
    return

file_num = 0
base_name = "N2_"

for i in range(101):
    dist = 0.8 + i*3.0/100.0
    file_num = file_num + 1
    WriteXYZFile(base_name+str(file_num-1)+".xyz",dist)
