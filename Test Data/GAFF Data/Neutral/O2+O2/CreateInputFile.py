import copy
import numpy as np

def WriteXYZFile(file_name,dist,angle):
    # Write the XYZ data to file_name based on the specified O...O distance and
    # Euler (z-axis) angle of the O2 dimer.
    
    OO_len = 1.16
    
    # Unrotated O2 Coordinates
    O1X = -OO_len
    O1Y = 0.0
    O1Z = 0.0
    
    O2X = 0.0
    O2Y = 0.0
    O2Z = 0.0
    
    O3X = dist
    O3Y = 0.0
    O3Z = 0.0
    
    O4X = OO_len + dist
    O4Y = 0.0
    O4Z = 0.0
    
    # Rotate the second O2 molecule
    O3X_rot = np.cos(angle)*O3X - np.sin(angle)*O3Y
    O3Y_rot = np.sin(angle)*O3X + np.cos(angle)*O3Y
    
    O4X_rot = np.cos(angle)*O4X - np.sin(angle)*O4Y
    O4Y_rot = np.sin(angle)*O4X + np.cos(angle)*O4Y
    
    fileID = open(file_name,"w+")
    fileID.write("4\n\n")
    
    fileID.write("O  ")
    fileID.write("{:12.6f} ".format(O1X))
    fileID.write("{:12.6f} ".format(O1Y))
    fileID.write("{:12.6f} ".format(O1Z)+"\n")
    
    fileID.write("O  ")
    fileID.write("{:12.6f} ".format(O2X))
    fileID.write("{:12.6f} ".format(O2Y))
    fileID.write("{:12.6f} ".format(O2Z)+"\n")
    
    fileID.write("O  ")
    fileID.write("{:12.6f} ".format(O3X_rot))
    fileID.write("{:12.6f} ".format(O3Y_rot))
    fileID.write("{:12.6f} ".format(O3Z)+"\n")
    
    fileID.write("O  ")
    fileID.write("{:12.6f} ".format(O4X_rot))
    fileID.write("{:12.6f} ".format(O4Y_rot))
    fileID.write("{:12.6f} ".format(O4Z)+"\n")
    
    fileID.write("\n\n")
    fileID.close()
    return

file_num = 0
base_name = "O2_O2_"

for i in range(101):
    dist = 1.5 + i*4.0/100.0
    for j in range(3):
        file_num = file_num + 1
        angle = j*(np.pi/(4.0))
        WriteXYZFile(base_name+str(file_num-1)+".xyz",dist,angle)
