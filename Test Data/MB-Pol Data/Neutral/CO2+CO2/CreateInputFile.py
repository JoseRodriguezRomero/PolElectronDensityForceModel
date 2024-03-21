import copy
import numpy as np

def WriteXYZFile(file_name,dist,angle):
    # Write the XYZ data to file_name based on the specified O...O distance and
    # Euler (z-axis) angle of the O2 dimer.
    
    CO_len = 1.16
    
    # Unrotated O2 Coordinates
    O1X = -2*CO_len
    O1Y = 0.0
    O1Z = 0.0
    
    C2X = -1*CO_len
    C2Y = 0.0
    C2Z = 0.0
    
    O3X = 0.0
    O3Y = 0.0
    O3Z = 0.0
    
    O4X = dist
    O4Y = 0.0
    O4Z = 0.0
    
    C5X = dist + 1*CO_len
    C5Y = 0.0
    C5Z = 0.0
    
    O6X = dist + 2*CO_len
    O6Y = 0.0
    O6Z = 0.0
    
    # Rotate the second CO2 molecule
    O4X_rot = np.cos(angle)*O4X - np.sin(angle)*O4Y
    O4Y_rot = np.sin(angle)*O4X + np.cos(angle)*O4Y
    
    C5X_rot = np.cos(angle)*C5X - np.sin(angle)*C5Y
    C5Y_rot = np.sin(angle)*C5X + np.cos(angle)*C5Y
    
    O6X_rot = np.cos(angle)*O6X - np.sin(angle)*O6Y
    O6Y_rot = np.sin(angle)*O6X + np.cos(angle)*O6Y
    
    fileID = open(file_name,"w+")
    fileID.write("SYSTEM CO2 CO2\n")
    fileID.write("MOLECULE\n")
    fileID.write("MONOMER CO2 \n")
    fileID.write("C  ")
    fileID.write("{:12.6f} ".format(C2X))
    fileID.write("{:12.6f} ".format(C2Y))
    fileID.write("{:12.6f} ".format(C2Z)+"\n")
    
    fileID.write("O  ")
    fileID.write("{:12.6f} ".format(O1X))
    fileID.write("{:12.6f} ".format(O1Y))
    fileID.write("{:12.6f} ".format(O1Z)+"\n")
    
    fileID.write("O  ")
    fileID.write("{:12.6f} ".format(O3X))
    fileID.write("{:12.6f} ".format(O3Y))
    fileID.write("{:12.6f} ".format(O3Z)+"\n")
    fileID.write("ENDMON \n")
    fileID.write("ENDMOL \n")
    
    fileID.write("MOLECULE\n")
    fileID.write("MONOMER CO2 \n")
    fileID.write("C  ")
    fileID.write("{:12.6f} ".format(C5X_rot))
    fileID.write("{:12.6f} ".format(C5Y_rot))
    fileID.write("{:12.6f} ".format(C5Z)+"\n")
    
    fileID.write("O  ")
    fileID.write("{:12.6f} ".format(O4X_rot))
    fileID.write("{:12.6f} ".format(O4Y_rot))
    fileID.write("{:12.6f} ".format(O4Z)+"\n")
    
    fileID.write("O  ")
    fileID.write("{:12.6f} ".format(O6X_rot))
    fileID.write("{:12.6f} ".format(O6Y_rot))
    fileID.write("{:12.6f} ".format(O6Z)+"\n")
    
    fileID.write("ENDMON \n")
    fileID.write("ENDMOL \n")
    
    fileID.write("ENDSYS \n\n")
    fileID.close()
    return

file_num = 0
base_name = "CO2_CO2_"

for i in range(101):
    dist = 1.5 + i*4.0/100.0
    for j in range(3):
        file_num = file_num + 1
        angle = j*(np.pi/(4.0))
        WriteXYZFile(base_name+str(file_num-1)+".nrg",dist,angle)
