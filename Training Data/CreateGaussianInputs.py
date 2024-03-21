import copy
import numpy as np

def RotateAndTranslate(coords,angs_and_disp):
    # Rotates and translates the coords by the angles and displacements
    # specified by angs_and_disp
    thx_x = float(angs_and_disp[0])
    thx_y = float(angs_and_disp[1])
    thx_z = float(angs_and_disp[2])
    dx = float(angs_and_disp[3])
    dy = float(angs_and_disp[4])
    dz = float(angs_and_disp[5])
    
    for i in range(len(coords)):
        old_coords = copy.copy(coords[i])
        coords[i][0] = np.cos(thx_z)*old_coords[0] - np.sin(thx_z)*old_coords[1]
        coords[i][1] = np.sin(thx_z)*old_coords[0] + np.cos(thx_z)*old_coords[1]
        
        old_coords = copy.copy(coords[i])
        coords[i][0] = np.cos(thx_y)*old_coords[0] + np.sin(thx_y)*old_coords[2]
        coords[i][2] = np.cos(thx_y)*old_coords[2] - np.sin(thx_y)*old_coords[0]

        old_coords = copy.copy(coords[i])
        coords[i][1] = np.cos(thx_x)*old_coords[1] - np.sin(thx_x)*old_coords[2]
        coords[i][2] = np.sin(thx_x)*old_coords[1] + np.cos(thx_x)*old_coords[2]
    
        coords[i][0] = coords[i][0] + dx*0.529177210903
        coords[i][1] = coords[i][1] + dy*0.529177210903
        coords[i][2] = coords[i][2] + dz*0.529177210903
    
    return coords
    
    
def ReadBaseCoords(file_name):
    # Reads the atom coordinates of the unrotated and untranslated molecule.
    coords = []
    fileID = open(file_name,"r")
    fileID.readline()
    
    while True:
        line = fileID.readline()
        line_splitted = line.split()
        
        if len(line_splitted) < 4:
            break
            
        new_coord = []
        for i in range(3):
            new_coord.append(0.529177210903*float(line_splitted[i]))
            
        coords.append(new_coord)
    return coords


def WriteGaussianFile(file_name,angs_and_disp):
    # Writes a formatted submission for Gaussian GD3-B3LYP/CEP-31G
    fileID = open(file_name + ".com","w")
    fileID.write("%NProcShared=4 \n")
    fileID.write("%Mem=4GB  \n")
    fileID.write("%Chk=" + file_name + ".chk \n\n")

    fileID.write("#p B3LYP/CEP-31G EmpiricalDispersion=GD3 \n\n")

    fileID.write(file_name + "\n\n")
    fileID.write("0 1\n")

    molec_a = ReadBaseCoords("H2O_ecp_fitted_data.txt")
    molec_b = ReadBaseCoords("H2O_ecp_fitted_data.txt")
    molec_b = RotateAndTranslate(molec_b,angs_and_disp)

    format_1 = "{:5s}";
    format_2 = "{:18.6f}";
    
    for i in range(2):
        if i == 0:
            molec = molec_a
        else:
            molec = molec_b
    
        fileID.write(format_1.format("O"))
        fileID.write(format_2.format(molec[0][0]))
        fileID.write(format_2.format(molec[0][1]))
        fileID.write(format_2.format(molec[0][2]))
        fileID.write("\n")
        
        fileID.write(format_1.format("H"))
        fileID.write(format_2.format(molec[1][0]))
        fileID.write(format_2.format(molec[1][1]))
        fileID.write(format_2.format(molec[1][2]))
        fileID.write("\n")
        
        fileID.write(format_1.format("H"))
        fileID.write(format_2.format(molec[2][0]))
        fileID.write(format_2.format(molec[2][1]))
        fileID.write(format_2.format(molec[2][2]))
        fileID.write("\n")
    
    fileID.write("\n")
    fileID.close()
    return


fileID = open("RotationsAndDisplacements.txt")
fileID.readline()

lines = fileID.readlines()
fileID.close()

base_name = "H2O_H2O_"

for i in range(len(lines)):
    line = lines[i]
    line_splitted = line.split()
    if len(line_splitted) > 1:
        file_name = base_name + str(i)
        WriteGaussianFile(file_name,line_splitted);
