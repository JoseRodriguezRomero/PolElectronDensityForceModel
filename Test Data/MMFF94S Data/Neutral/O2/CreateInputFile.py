import copy
import numpy as np

def WriteXYZFile(file_name,dist):
    # Write the XYZ data to file_name based on the specified O...O interatomic
    # O-O bond distance.
    
    O1X = 0.0;
    O1Y = 0.0;
    O1Z = 0.0;
    
    O2X = dist;
    O2Y = 0.0;
    O2Z = 0.0;
    
    fileID = open(file_name,"w+")
    fileID.write("2\n\n")
    
    fileID.write("O  ")
    fileID.write("{:12.6f} ".format(O1X))
    fileID.write("{:12.6f} ".format(O1Y))
    fileID.write("{:12.6f} ".format(O1Z)+"\n")
    
    fileID.write("O  ")
    fileID.write("{:12.6f} ".format(O2X))
    fileID.write("{:12.6f} ".format(O2Y))
    fileID.write("{:12.6f} ".format(O2Z)+"\n")
    
    fileID.write("\n\n")
    fileID.close()
    return

file_num = 0
base_name = "O2_"

for i in range(101):
    dist = 0.8 + i*3.0/100.0
    file_num = file_num + 1
    WriteXYZFile(base_name+str(file_num-1)+".xyz",dist)
