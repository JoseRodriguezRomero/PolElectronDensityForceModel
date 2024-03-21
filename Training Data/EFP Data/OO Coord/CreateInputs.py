def WriteGaussianFile(file_name,angs_and_disp):
    # Writes a formatted submission for Gaussian GD3-B3LYP/CEP-31G
    fileID = open(file_name + ".inp","w")
    fileID.write("fragment h2o_l \n")
    fileID.write("   ");
    
    for i in range(6):
        fileID.write("{:18.6f}".format(0.0) + " ")

    fileID.write("\n")
    fileID.write("fragment h2o_l \n")
    fileID.write("   ");
    for i in range(3):
        aux_num = float(angs_and_disp[3+i])*0.529177210903
        fileID.write("{:18.6f}".format(aux_num) + " ")
    for i in range(3):
        aux_num = float(angs_and_disp[i])*0.529177210903
        fileID.write("{:18.6f}".format(aux_num) + " ")
        
    fileID.write("\n")
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
