import numpy as np

def WriteInputFile(file_name,dist,angle):
    # Writes a formatted submission for Gaussian GD3-B3LYP/CEP-31G
    fileID = open(file_name + ".inp","w")
    fileID.write("fragment h2_l \n")
    fileID.write("   ");
    
    HH_len = 0.74
    
    fileID.write("{:18.6f}".format(-HH_len/2) + " ")
    fileID.write("{:18.6f}".format(0.0) + " ")
    fileID.write("{:18.6f}".format(0.0) + " ")
    fileID.write("{:18.6f}".format(0.0) + " ")
    fileID.write("{:18.6f}".format(0.0) + " ")
    fileID.write("{:18.6f}".format(0.0) + " ")
        
    fileID.write("\n")
    fileID.write("fragment h2_l \n")
    fileID.write("   ");
    fileID.write("{:18.6f}".format((dist+HH_len/2)*np.cos(angle)) + " ")
    fileID.write("{:18.6f}".format((dist+HH_len/2)*np.sin(angle)) + " ")
    fileID.write("{:18.6f}".format(0.0) + " ")
    fileID.write("{:18.6f}".format(0.0) + " ")
    fileID.write("{:18.6f}".format(0.0) + " ")
    fileID.write("{:18.6f}".format(angle) + " ")
            
    fileID.write("\n")
    return

file_num = 0
base_name = "H2_H2_"

for i in range(101):
    dist = 1.0 + i*4.0/100.0
    for j in range(3):
        file_num = file_num + 1
        angle = j*(np.pi/(4.0))
        WriteInputFile(base_name+str(file_num-1),dist,angle)
