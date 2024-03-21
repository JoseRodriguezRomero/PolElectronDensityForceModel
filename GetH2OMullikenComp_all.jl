using Plots, SpecialFunctions, LaTeXStrings, StatsPlots, StatsBase, Statistics, PGFPlotsX
include("GetPotentialFromAngles.jl")

function GetEnergyFromLogFile(file_name::String)
    # Reads the energy (in Hartree) of the molecule from the log file specified
    # by file_name.
    energy = 0.0;
    fileID = open(file_name,"r");
 
    for line in readlines(fileID)
        if contains(line,"E(RB3LYP)") || contains(line,"E(UB3LYP)")
            line_splitted = split(line);
            energy = parse(Float64,line_splitted[5]);
        end

        if contains(line,"CCSD(T)= ")
            aux_str = split(line)[end];
            energy = parse(Float64,replace(aux_str, "D" => "E", count = 1));
        end
    end

    close(fileID);
    return energy;
end

function GetGeometryFromLogFile(file_name::String)
    # Reads the geometry of the molecule from the Gaussian logfile specified by 
    # file_name.
    geometry = zeros(Float64,0,0);
    fileID = open(file_name,"r");

    lines = readlines(fileID);
    for i in 1:lastindex(lines)
        line = lines[i];

        if contains(line,"Standard orientation:")
            j = i + 5;
            geometry = zeros(Float64,0,3);
            while true
                line = lines[j];
                if contains(line,"----------------------------------")
                    break;
                end

                line_splitted = split(line);
                x = parse(Float64,line_splitted[4]);
                y = parse(Float64,line_splitted[5]);
                z = parse(Float64,line_splitted[6]);
                geometry = [geometry; [x,y,z]'];

                j += 1;
            end
        end
    end

    close(fileID);
    return geometry;
end

function ReadGaussianMullikenCharges(file_name::String)
    # Returns a vector with the Mulliken charges of each atom, in the 
    # same order as given in the Gaussian input file, from the Gaussian 
    # logfile.
    charges = zeros(Float64,0);

    fileID = open(file_name);
    lines = readlines(fileID);

    for i in eachindex(lines)
        line = lines[i];
        
        if contains(line,"Mulliken charges:")
            j = 2;
            while true
                line = lines[i+j];

                if contains(line,"Sum of Mulliken charges =")
                    break
                end

                line_splitted = split(line);
                charge =  parse(Float64,line_splitted[end]);
                push!(charges,charge)
                j += 1;
            end

            break; 
        end
    end

    close(fileID);
    return charges;
end

function ReadXCCoeffs(file_name::String)
    # Reads the XC coefficients from the specified text file. The file is 
    # assumed to be in a single column format.
    ret = zeros(Float64,0);
    for line in readlines(file_name)
        push!(ret,parse(Float64,line));
    end

    return ret;
end

ecp_mull_charges = zeros(Float64,10000,6);
model_ecp_mull_charges = zeros(Float64,10000,6);

fulle_mull_charges = zeros(Float64,10000,6);
model_fulle_mull_charges = zeros(Float64,10000,6);

model1_xc_coeffs_ECP = ReadXCCoeffs("ECP_XCCoeffs.txt");
model1_xc_coeffs_FullE = ReadXCCoeffs("FullE_XCCoeffs.txt");
model1_xc_coeffs_ECP_Pol = ReadXCCoeffs("ECP_XCCoeffs_Pol.txt");
model1_xc_coeffs_FullE_Pol = ReadXCCoeffs("FullE_XCCoeffs_Pol.txt");

ECP_order = ceil(Int,length(model1_xc_coeffs_ECP)/4.0) - 1;
FullE_order = ceil(Int,length(model1_xc_coeffs_FullE)/4.0) - 1;

bond_θ = -(π/2.0)*(13.0/25.0);

OH_angs_and_disps = zeros(Float64,5000,6);
OO_angs_and_disps = zeros(Float64,5000,6);

lines = readlines("Training Data/Gaussian Data/FullE/OH Coord/RotationsAndDisplacements.txt");
for i in eachindex(lines[2:end])
    line_splitted = split(lines[i+1]);
    for j in eachindex(line_splitted)
        OH_angs_and_disps[i,j] = parse(Float64,line_splitted[j]);
    end
end

lines = readlines("Training Data/Gaussian Data/FullE/OO Coord/RotationsAndDisplacements.txt");
for i in eachindex(lines[2:end])
    line_splitted = split(lines[i+1]);
    for j in eachindex(line_splitted)
        OO_angs_and_disps[i,j] = parse(Float64,line_splitted[j]);
    end
end

ecp_mol_a = ReadMolecule("H2O_ecp_fitted_data.txt");
fulle_mol_a = ReadMolecule("H2O_fullE_fitted_data.txt");

ecp_h2o_energy = GetEnergyFromLogFile("Training Data/Gaussian Data/ECP/H2O.log");
fulle_h2o_energy = GetEnergyFromLogFile("Training Data/Gaussian Data/FullE/H2O.log");

ecp_energies = zeros(Float64,10000);
fulle_energies = zeros(Float64,10000);
for i in 0:4999
    base_name = "Training Data/Gaussian Data/ECP/OH Coord/H2O_H2O_";
    local file_name = base_name * string(i) * ".log";
    ecp_energies[i+1] = GetEnergyFromLogFile(file_name) - 2.0*ecp_h2o_energy;

    base_name = "Training Data/Gaussian Data/ECP/OO Coord/H2O_H2O_";
    local file_name = base_name * string(i) * ".log";
    ecp_energies[i+5001] = GetEnergyFromLogFile(file_name) - 2.0*ecp_h2o_energy;

    base_name = "Training Data/Gaussian Data/FullE/OH Coord/H2O_H2O_";
    local file_name = base_name * string(i) * ".log";
    fulle_energies[i+1] = GetEnergyFromLogFile(file_name) - 2.0*fulle_h2o_energy;

    base_name = "Training Data/Gaussian Data/FullE/OO Coord/H2O_H2O_";
    local file_name = base_name * string(i) * ".log";
    fulle_energies[i+5001] = GetEnergyFromLogFile(file_name) - 2.0*fulle_h2o_energy;
end

kj_per_mol = 2625.5002;
ecp_energies .*= kj_per_mol;
fulle_energies .*= kj_per_mol;

for i in 0:4999
    # OH Coord
    base_name = "Training Data/Gaussian Data/FullE/OH Coord/H2O_H2O_";
    local file_name = base_name * string(i) * ".log";
    new_charges = ReadGaussianMullikenCharges(file_name);

    mol_b = copy(fulle_mol_a);
    MoveAndRotateMolec!(mol_b,OH_angs_and_disps[i+1,:]);
    PolarizeMolecules!([fulle_mol_a,mol_b],model1_xc_coeffs_FullE_Pol);

    num_atoms = size(fulle_mol_a.atoms_data)[1];
    num_clouds = size(fulle_mol_a.cloud_data)[1];
    clouds_per_atom = ceil(Int,num_clouds/num_atoms);

    fulle_mull_charges[i+1,:] = new_charges;
    
    ρa = fulle_mol_a.atoms_data[:,4];
    ρb = mol_b.atoms_data[:,4];

    for i in 1:clouds_per_atom
        new_ρ = fulle_mol_a.cloud_data[i:clouds_per_atom:end,4];
        new_ρ .*= fulle_mol_a.cloud_data[i:clouds_per_atom:end,6];
        ρa -= new_ρ;

        new_ρ = mol_b.cloud_data[i:clouds_per_atom:end,4];
        new_ρ .*= mol_b.cloud_data[i:clouds_per_atom:end,6];
        ρb -= new_ρ;
    end

    model_fulle_mull_charges[i+1,:] = vcat(ρa,ρb);

    # OO Coord
    base_name = "Training Data/Gaussian Data/FullE/OO Coord/H2O_H2O_";
    local file_name = base_name * string(i) * ".log" 
    mol_geom = GetGeometryFromLogFile(file_name);
    new_charges = ReadGaussianMullikenCharges(file_name);

    mol_b = copy(fulle_mol_a);
    MoveAndRotateMolec!(mol_b,OO_angs_and_disps[i+1,:]);
    PolarizeMolecules!([fulle_mol_a,mol_b],model1_xc_coeffs_FullE_Pol);

    num_atoms = size(fulle_mol_a.atoms_data)[1];
    num_clouds = size(fulle_mol_a.cloud_data)[1];
    clouds_per_atom = ceil(Int,num_clouds/num_atoms);

    fulle_mull_charges[i+5001,:] = new_charges;

    ρa = fulle_mol_a.atoms_data[:,4];
    ρb = mol_b.atoms_data[:,4];

    for i in 1:clouds_per_atom
        new_ρ = fulle_mol_a.cloud_data[i:clouds_per_atom:end,4];
        new_ρ .*= fulle_mol_a.cloud_data[i:clouds_per_atom:end,6];
        ρa -= new_ρ;

        new_ρ = mol_b.cloud_data[i:clouds_per_atom:end,4];
        new_ρ .*= mol_b.cloud_data[i:clouds_per_atom:end,6];
        ρb -= new_ρ;
    end

    model_fulle_mull_charges[i+5001,:] = vcat(ρa,ρb);
end

for i in 0:4999
    # OH Coord
    base_name = "Training Data/Gaussian Data/ECP/OH Coord/H2O_H2O_";
    local file_name = base_name * string(i) * ".log";
    new_charges = ReadGaussianMullikenCharges(file_name);

    mol_b = copy(ecp_mol_a);
    MoveAndRotateMolec!(mol_b,OH_angs_and_disps[i+1,:]);
    PolarizeMolecules!([ecp_mol_a,mol_b],model1_xc_coeffs_ECP_Pol);

    num_atoms = size(ecp_mol_a.atoms_data)[1];
    num_clouds = size(ecp_mol_a.cloud_data)[1];
    clouds_per_atom = ceil(Int,num_clouds/num_atoms);

    ecp_mull_charges[i+1,:] = new_charges;
    
    ρa = ecp_mol_a.atoms_data[:,4];
    ρb = mol_b.atoms_data[:,4];

    for i in 1:clouds_per_atom
        new_ρ = ecp_mol_a.cloud_data[i:clouds_per_atom:end,4];
        new_ρ .*= ecp_mol_a.cloud_data[i:clouds_per_atom:end,6];
        ρa -= new_ρ;

        new_ρ = mol_b.cloud_data[i:clouds_per_atom:end,4];
        new_ρ .*= mol_b.cloud_data[i:clouds_per_atom:end,6];
        ρb -= new_ρ;
    end

    model_ecp_mull_charges[i+1,:] = vcat(ρa,ρb);

    # OO Coord
    base_name = "Training Data/Gaussian Data/ECP/OO Coord/H2O_H2O_";
    local file_name = base_name * string(i) * ".log" 
    mol_geom = GetGeometryFromLogFile(file_name);

    mol_b = copy(ecp_mol_a);
    MoveAndRotateMolec!(mol_b,OO_angs_and_disps[i+1,:]);
    PolarizeMolecules!([ecp_mol_a,mol_b],model1_xc_coeffs_ECP_Pol);

    num_atoms = size(ecp_mol_a.atoms_data)[1];
    num_clouds = size(ecp_mol_a.cloud_data)[1];
    clouds_per_atom = ceil(Int,num_clouds/num_atoms);

    ecp_mull_charges[i+5001,:] = new_charges;
    
    ρa = ecp_mol_a.atoms_data[:,4];
    ρb = mol_b.atoms_data[:,4];

    for i in 1:clouds_per_atom
        new_ρ = ecp_mol_a.cloud_data[i:clouds_per_atom:end,4];
        new_ρ .*= ecp_mol_a.cloud_data[i:clouds_per_atom:end,6];
        ρa -= new_ρ;

        new_ρ = mol_b.cloud_data[i:clouds_per_atom:end,4];
        new_ρ .*= mol_b.cloud_data[i:clouds_per_atom:end,6];
        ρb -= new_ρ;
    end

    model_ecp_mull_charges[i+5001,:] = vcat(ρa,ρb);
end

# ECP Plot
p1 = plot([-2, 2],[-2, 2]);
scatter!(ecp_mull_charges[ecp_energies .>= 50.0,:][:], 
    model_ecp_mull_charges[ecp_energies .>= 50.0,:][:],
    markersize=1.5,legend=false, markerstrokewidth=0.05);
scatter!(ecp_mull_charges[ecp_energies .< 50.0,:][:], 
    model_ecp_mull_charges[ecp_energies .< 50.0,:][:],
    markersize=1.5,legend=false, markerstrokewidth=0.2);
plot!(xlims=(-1,0.6),ylims=(-1,0.6));
plot!(xticks=-1:0.4:0.6,yticks=-1:0.4:0.6);
plot!(xlabel="Mulliken Charges\nCCSD(T)/CEP-31G");
plot!(ylabel="Modeled\nPartial Charges");
plot!(xguidefontsize=8);
plot!(yguidefontsize=8);

ss_res = sum((ecp_mull_charges[:] - model_ecp_mull_charges[:]).^2.0);
ss_tot = sum((ecp_mull_charges[:] .- mean(ecp_mull_charges[:])).^2.0);
R² = 1.0 - ss_res/ss_tot;
ecp_R² = 1.0 - ss_res/ss_tot;

aux_text = "R² = "*string(round(R²*10000000)/10000000);
annotate!(0.55,-0.885,text(aux_text,:center,:right,8));

# Full Electron Plot
p2 = plot([-2, 2],[-2, 2]);
scatter!(fulle_mull_charges[fulle_energies .>= 50.0,:][:], 
    model_fulle_mull_charges[fulle_energies .>= 50.0,:][:],
    markersize=1.5,legend=false, markerstrokewidth=0.05);
scatter!(fulle_mull_charges[fulle_energies .< 50.0,:][:], 
    model_fulle_mull_charges[fulle_energies .< 50.0,:][:],
    markersize=1.5,legend=false, markerstrokewidth=0.2);
plot!(xlims=(-1,0.6),ylims=(-1,0.6));
plot!(xticks=-1:0.4:0.6,yticks=(-1:0.4:0.6,[]));
plot!(ylabel="");
plot!(xlabel="Mulliken Charges\nCCSD(T)/cc-pVTZ");
plot!(xguidefontsize=8);
plot!(yguidefontsize=8);

ss_res = sum((fulle_mull_charges[:] - model_fulle_mull_charges[:]).^2.0);
ss_tot = sum((fulle_mull_charges[:] .- mean(fulle_mull_charges[:])).^2.0);
R² = 1.0 - ss_res/ss_tot;
fulle_R² = 1.0 - ss_res/ss_tot;

aux_text = "R² = "*string(round(R²*10000000)/10000000);
annotate!(0.55,-0.885,text(aux_text,:center,:right,8));

plot(p1,p2);
plot!(left_margin=3Plots.mm);
plot!(right_margin=2Plots.mm);
plot!(bottom_margin=4Plots.mm);
plot!(top_margin=1Plots.mm);
plot!(size=(650,240))

plot!(dpi=1000);
savefig("mulliken_comps.png");

PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.DEFAULT_PREAMBLE; "\\usepackage{amsmath,accents}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage{accents}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage[utf8]{inputenc}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage[T1]{fontenc}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage{fontspec}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage{mathptmx}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage[fontsize=10pt]{fontsize}"];

# PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.DEFAULT_PREAMBLE; "\\usepackage[bitstream-charter]{mathdesign}"];
# PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage{amsmath,mathtools,accents}"];
# PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage{accents}"];
# PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage[utf8]{inputenc}"];
# PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage[T1]{fontenc}"];
# PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage{fontspec}"];
# PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage[fontsize=10pt]{fontsize}"];

fulle_low_x = fulle_mull_charges[fulle_energies .< 50.0,:][:];
fulle_low_y = model_fulle_mull_charges[fulle_energies .< 50.0,:][:];
fulle_high_x = fulle_mull_charges[fulle_energies .>= 50.0,:][:];
fulle_high_y = model_fulle_mull_charges[fulle_energies .>= 50.0,:][:];

ecp_low_x = ecp_mull_charges[ecp_energies .< 50.0,:][:];
ecp_low_y = model_ecp_mull_charges[ecp_energies .< 50.0,:][:];
ecp_high_x = ecp_mull_charges[ecp_energies .>= 50.0,:][:];
ecp_high_y = model_ecp_mull_charges[ecp_energies .>= 50.0,:][:];

dummy_inds = Int.(round.(60000.0.*rand(30)));
dummy_x = fulle_mull_charges[dummy_inds];
dummy_y = model_fulle_mull_charges[dummy_inds];

hist_bins = collect(-1.0:0.001:0.6);
hist_ecp = zeros(Float64,length(hist_bins));
hist_model_ecp = zeros(Float64,length(hist_bins));
hist_fulle = zeros(Float64,length(hist_bins));
hist_model_fulle = zeros(Float64,length(hist_bins));

for val in ecp_mull_charges[:]
    bin = Int(floor((val + 1.0)/0.001));
    hist_ecp[bin] = hist_ecp[bin] + 1;
end

for val in model_ecp_mull_charges[:]
    bin = Int(floor((val + 1.0)/0.001));
    hist_model_ecp[bin] = hist_model_ecp[bin] + 1;
end

for val in fulle_mull_charges[:]
    bin = Int(floor((val + 1.0)/0.001));
    hist_fulle[bin] = hist_fulle[bin] + 1;
end

for val in model_fulle_mull_charges[:]
    bin = Int(floor((val + 1.0)/0.001));
    hist_model_fulle[bin] = hist_model_fulle[bin] + 1;
end

p1 = @pgf Axis(
    {
        xmajorgrids,
        ymajorgrids,
        width = 175, 
        height = 140,
        xmin = -1,
        xmax = 0.6,
        ymin = -1,
        ymax = 0.6,
        ytick = "-1,-0.6,-0.2,0.2,0.6",
        xtick = "-1,-0.6,-0.2,0.2,0.6",
        xticklabels = L"$-1.0$,$-0.6$,$-0.2$,$0.2$,$0.6$",
        yticklabels = L"$-1.0$,$-0.6$,$-0.2$,$0.2$,$0.6$",
        ylabel = L"\begin{gathered} \mathrm{Modeled} \\[-0.15cm] \mathrm{Partial \ Charges} \end{gathered}",
        xlabel = L"\begin{gathered} \mathrm{Mulliken \ Charges} \\[-0.15cm] \mathrm{CCSD(T)/cc\text{-}pVTZ} \end{gathered}",
        "grid style={line width=.1pt, draw=gray!10},major grid style={line width=.2pt,draw=gray!25}",
    },
    Plot(
        {
            "only marks",
            mark = "x",
            "mark options" = {"color" = theme_palette(:auto)[2], scale = 0.2},
        },
        Table([fulle_high_x[:],fulle_high_y[:]]),
        "node[] at (-0.6,0.45) {\\fontsize{8pt}{8pt}\\selectfont \$R^2 = "*string(round(fulle_R²,digits=5))*"\$}"
    ),
    Plot(
        {
            "only marks",
            mark = "x",
            "mark options" = {"color" = theme_palette(:auto)[3], scale = 0.2},
        },
        Table([fulle_low_x[:],fulle_low_y[:]])
    ),
    Plot(
        {
            no_markers,
            color = theme_palette(:auto)[1]
        },
        Coordinates([-2,2],[-2,2])
    ),
);

pg11 = @pgf Axis(
    {
        width = 175, 
        height = 80,
        xmin = -1,
        xmax = 0.6,
        ymin = 0,
        ymax = 0.15,
        yticklabels = "\\empty",
        xticklabels = "\\empty",
        "hide axis"
    },
    Plot(
        {
            no_marks,
        },
        Coordinates(hist_bins,hist_fulle ./ 60000.0),
    ),
);

pg12 = @pgf Axis(
    {
        width = 140, 
        height = 140,
        ymin = -1,
        ymax = 0.6,
        xmin = 0,
        xmax = 0.15*(140.0/80.0),
        yticklabels = "\\empty",
        xticklabels = "\\empty",
        "hide axis"
    },
    Plot(
        {
            no_marks,
        },
        Coordinates(hist_model_fulle ./ 60000.0, hist_bins),
    ),
);

pg10 = @pgf Axis(
    {
        height = 80, 
        width = 140, 
        ymin = -1,
        ymax = 0.6,
        xmin = 0,
        xmax = 0.15*(140.0/80.0),
        yticklabels = "\\empty",
        xticklabels = "\\empty",
        "hide axis"
    },
    Plot(
        {
            no_marks,
        },
        Coordinates([0],[0]),
    ),
);

# @pgf GroupPlot({group_style = { group_size = "2 by 2","vertical sep = 0","horizontal sep = 0"}}, pg11, pg0, p1, pg12)

p2 = @pgf Axis(
    {
        xmajorgrids,
        ymajorgrids,
        width = 175, 
        height = 140,
        xmin = -1,
        xmax = 0.6,
        ymin = -1,
        ymax = 0.6,
        ytick = "-1,-0.6,-0.2,0.2,0.6",
        xtick = "-1,-0.6,-0.2,0.2,0.6",
        xticklabels = L"$-1.0$,$-0.6$,$-0.2$,$0.2$,$0.6$",
        yticklabels = L"$-1.0$,$-0.6$,$-0.2$,$0.2$,$0.6$",
        yticklabels = "\\empty",
        xlabel = L"\begin{gathered} \mathrm{Mulliken \ Charges} \\[-0.15cm] \mathrm{CCSD(T)/CEP\text{-}31G} \end{gathered}",
        "grid style={line width=.1pt, draw=gray!10},major grid style={line width=.2pt,draw=gray!25}",
    },
    Plot(
        {
            "only marks",
            mark = "x",
            "mark options" = {"color" = theme_palette(:auto)[2], scale = 0.2},
        },
        Table([ecp_high_x[:],ecp_high_y[:]]),
        "node[] at (-0.6,0.45) {\\fontsize{8pt}{8pt}\\selectfont \$R^2 = "*string(round(ecp_R²,digits=5))*"\$}"
    ),
    Plot(
        {
            "only marks",
            mark = "x",
            "mark options" = {"color" = theme_palette(:auto)[3], scale = 0.2},
        },
        Table([ecp_low_x[:],ecp_low_y[:]])
    ),
    Plot(
        {
            no_markers,
            color = theme_palette(:auto)[1]
        },
        Coordinates([-2,2],[-2,2]),
    ),
);

pg21 = @pgf Axis(
    {
        width = 175, 
        height = 80,
        xmin = -1,
        xmax = 0.6,
        ymin = 0,
        ymax = 0.15,
        yticklabels = "\\empty",
        xticklabels = "\\empty",
        "hide axis"
    },
    Plot(
        {
            no_marks,
        },
        Coordinates(hist_bins,hist_ecp ./ 60000.0),
    ),
);

pg22 = @pgf Axis(
    {
        width = 80, 
        height = 140,
        ymin = -1,
        ymax = 0.6,
        xmin = 0,
        xmax = 0.15,
        yticklabels = "\\empty",
        xticklabels = "\\empty",
        "hide axis"
    },
    Plot(
        {
            no_marks,
        },
        Coordinates(hist_model_ecp ./ 60000.0, hist_bins),
    ),
);

pg20 = @pgf Axis(
    {
        height = 80, 
        width = 80, 
        ymin = -1,
        ymax = 0.6,
        xmin = 0,
        xmax = 0.2,
        yticklabels = "\\empty",
        xticklabels = "\\empty",
        "hide axis"
    },
    Plot(
        {
            no_marks,
        },
        Coordinates([0],[0]),
    ),
);

# @pgf GroupPlot({group_style = { group_size = "2 by 1",}}, p1, p2)
@pgf GroupPlot({group_style = { group_size = "4 by 4",
    "vertical sep = 0","horizontal sep = 0"}}, 
    pg11, pg10, pg21, pg20, p1, pg12, p2, pg22)

