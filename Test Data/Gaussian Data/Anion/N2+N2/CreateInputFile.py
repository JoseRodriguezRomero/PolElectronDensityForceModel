import copy
import numpy as np

def WriteXYZFile(file_name,dist,angle):
    # Write the XYZ data to file_name based on the specified O...O distance and
    # Euler (z-axis) angle of the O2 dimer.
    
    NN_len = 1.09
    
    # Unrotated O2 Coordinates
    N1X = -NN_len;
    N1Y =  0.0;
    N1Z =  0.0;
    
    N2X = 0.0;
    N2Y = 0.0;
    N2Z = 0.0;
    
    N3X = dist;
    N3Y = 0.0;
    N3Z = 0.0;
    
    N4X = NN_len + dist;
    N4Y = 0.0;
    N4Z = 0.0;
    
    # Rotate the second O2 molecule
    N3X_rot = np.cos(angle)*N3X - np.sin(angle)*N3Y
    N3Y_rot = np.sin(angle)*N3X + np.cos(angle)*N3Y
    
    N4X_rot = np.cos(angle)*N4X - np.sin(angle)*N4Y
    N4Y_rot = np.sin(angle)*N4X + np.cos(angle)*N4Y
    
    fileID = open(file_name+".com","w+")
    fileID.write("%NProcShared=4 \n")
    fileID.write("%Mem=4GB  \n")
    fileID.write("%Chk="+file_name+".chk \n\n")
    fileID.write("#p ROCCSD(T,MaxCyc=200)/cc-pVTZ SCF(MaxCyc=200)\n\n")
    fileID.write(file_name+" \n\n")
    fileID.write("-1 2\n")
    
    fileID.write("N  ")
    fileID.write("{:12.6f} ".format(N1X))
    fileID.write("{:12.6f} ".format(N1Y))
    fileID.write("{:12.6f} ".format(N1Z)+"\n")
    
    fileID.write("N  ")
    fileID.write("{:12.6f} ".format(N2X))
    fileID.write("{:12.6f} ".format(N2Y))
    fileID.write("{:12.6f} ".format(N2Z)+"\n")
    
    fileID.write("N  ")
    fileID.write("{:12.6f} ".format(N3X_rot))
    fileID.write("{:12.6f} ".format(N3Y_rot))
    fileID.write("{:12.6f} ".format(N3Z)+"\n")
    
    fileID.write("N  ")
    fileID.write("{:12.6f} ".format(N4X_rot))
    fileID.write("{:12.6f} ".format(N4Y_rot))
    fileID.write("{:12.6f} ".format(N4Z)+"\n")
    
    fileID.write("\n\n")
    fileID.close()
    return

file_num = 0
base_name = "N2_N2_"

for i in range(101):
    dist = 1.0 + i*4.0/100.0
    for j in range(3):
        file_num = file_num + 1
        angle = j*(np.pi/(4.0))
        WriteXYZFile(base_name+str(file_num-1),dist,angle)
