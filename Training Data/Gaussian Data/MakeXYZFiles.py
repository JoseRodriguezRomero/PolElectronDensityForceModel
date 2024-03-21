import numpy as np


def ElementFromZ(atomic_number):
    if atomic_number == 1:
        return "H"
    elif atomic_number == 2:
        return "He"
    elif atomic_number == 3:
        return "Li"
    elif atomic_number == 4:
        return "Be"
    elif atomic_number == 5:
        return "B"
    elif atomic_number == 6:
        return "C"
    elif atomic_number == 7:
        return "N"
    elif atomic_number == 8:
        return "O"
    elif atomic_number == 9:
        return "F"
    elif atomic_number == 10:
        return "Ne"
    elif atomic_number == 11:
        return "Na"
    elif atomic_number == 12:
        return "Mg"
    elif atomic_number == 13:
        return "Al"
    elif atomic_number == 14:
        return "Si"
    elif atomic_number == 15:
        return "P"
    elif atomic_number == 16:
        return "S"
    elif atomic_number == 17:
        return "Cl"
    elif atomic_number == 18:
        return "Ar"
    elif atomic_number == 19:
        return "K"
    elif atomic_number == 20:
        return "Ca"
    elif atomic_number == 21:
        return "Sc"
    elif atomic_number == 22:
        return "Ti"
    elif atomic_number == 23:
        return "V"
    elif atomic_number == 24:
        return "Cr"
    elif atomic_number == 25:
        return "Mn"
    elif atomic_number == 26:
        return "Fe"
    elif atomic_number == 27:
        return "Co"
    elif atomic_number == 28:
        return "Ni"
    elif atomic_number == 29:
        return "Cu"
    elif atomic_number == 30:
        return "Zn"
    
    # oh noes! I has the lazy!
    return "X"

def WriteFile(file_name):
    # Writes an XYZ file from the geometry of the Gaussian logfile.
    fileID1 = open(file_name+".log","r")
    fileID2 = open(file_name+".xyz","w")
    
    natoms = 0
    lines = fileID1.readlines()
    
    for i in range(len(lines)):
        line = lines[i]
        if "NAtoms=" in line:
            line_splitted = line.split()
            natoms = int(line_splitted[1])
            fileID2.write(str(natoms)+"\n\n")
        
        if "Input orientation:" in line:
            for j in range(natoms):
                line_splitted = lines[i+j+5].split()
                atomZ = int(line_splitted[1])
                atom_px = float(line_splitted[3])
                atom_py = float(line_splitted[4])
                atom_pz = float(line_splitted[5])
                fileID2.write("{:5}".format(ElementFromZ(atomZ)))
                fileID2.write("{:18.6f}".format(atom_px))
                fileID2.write("{:18.6f}".format(atom_py))
                fileID2.write("{:18.6f}".format(atom_pz))
                fileID2.write("\n")
            fileID2.write("\n")
            return
    
    fileID1.close()
    fileID2.close()
    return



base_dir = "OH Coord/"
WriteFile(base_dir + "H2O")
for j in range(5000):
    file_name = "H2O_H2O_" + str(j)
    WriteFile(base_dir + file_name)

base_dir = "OO Coord/"
WriteFile(base_dir + "H2O")
for j in range(5000):
    file_name = "H2O_H2O_" + str(j)
    WriteFile(base_dir + file_name)
