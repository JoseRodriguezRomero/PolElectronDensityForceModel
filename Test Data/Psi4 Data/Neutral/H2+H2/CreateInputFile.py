import copy
import numpy as np

def WriteXYZFile(file_name,dist,angle):
    # Write the XYZ data to file_name based on the specified O...O distance and
    # Euler (z-axis) angle of the O2 dimer.
    
    HH_len = 0.74
    
    # Unrotated O2 Coordinates
    H1X = -HH_len;
    H1Y =  0.0;
    H1Z =  0.0;
    
    H2X = 0.0;
    H2Y = 0.0;
    H2Z = 0.0;
    
    H3X = dist;
    H3Y = 0.0;
    H3Z = 0.0;
    
    H4X = HH_len + dist;
    H4Y = 0.0;
    H4Z = 0.0;
    
    # Rotate the second O2 molecule
    H3X_rot = np.cos(angle)*H3X - np.sin(angle)*H3Y
    H3Y_rot = np.sin(angle)*H3X + np.cos(angle)*H3Y
    
    H4X_rot = np.cos(angle)*H4X - np.sin(angle)*H4Y
    H4Y_rot = np.sin(angle)*H4X + np.cos(angle)*H4Y
    
    fileID = open(file_name+".com","w+")
    fileID.write("%NProcShared=4 \n")
    fileID.write("%Mem=4GB  \n")
    fileID.write("%Chk="+file_name+".chk \n\n")
    fileID.write("#p CCSD(T)/cc-pVTZ \n\n")
    fileID.write(file_name+" \n\n")
    fileID.write("0 1\n")
    
    fileID.write("H  ")
    fileID.write("{:12.6f} ".format(H1X))
    fileID.write("{:12.6f} ".format(H1Y))
    fileID.write("{:12.6f} ".format(H1Z)+"\n")
    
    fileID.write("H  ")
    fileID.write("{:12.6f} ".format(H2X))
    fileID.write("{:12.6f} ".format(H2Y))
    fileID.write("{:12.6f} ".format(H2Z)+"\n")
    
    fileID.write("H  ")
    fileID.write("{:12.6f} ".format(H3X_rot))
    fileID.write("{:12.6f} ".format(H3Y_rot))
    fileID.write("{:12.6f} ".format(H3Z)+"\n")
    
    fileID.write("H  ")
    fileID.write("{:12.6f} ".format(H4X_rot))
    fileID.write("{:12.6f} ".format(H4Y_rot))
    fileID.write("{:12.6f} ".format(H4Z)+"\n")
    
    fileID.write("\n\n")
    fileID.close()
    return

file_num = 0
base_name = "H2_H2_"

for i in range(101):
    dist = 1.0 + i*4.0/100.0
    for j in range(3):
        file_num = file_num + 1
        angle = j*(np.pi/(4.0))
        WriteXYZFile(base_name+str(file_num-1),dist,angle)
