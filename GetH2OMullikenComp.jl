using Plots, SpecialFunctions, LaTeXStrings, StatsPlots, PGFPlotsX
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

function ChargesAtDistanceOO(dist::Float64,pol_xc_coeffs::Vector{Float64},
    density_file::String)
    # Returns the energy of the uncorelated model of two molecules of water
    # on the OH reaction coordinate at distance dist.
    # order = trunc(Int,length(xc_coeffs)/4) - 1;
    a = ReadMolecule(density_file,Float64);
    b = copy(a);

    # CenterAtAtomIndex!(a,1);
    # CenterAtAtomIndex!(b,1);

    disp_and_angs = zeros(Float64,6);
    disp_and_angs[3] = π;
    disp_and_angs[5] = dist;

    MoveAndRotateMolec!(a,disp_and_angs);
    PolarizeMolecules!([a,b],pol_xc_coeffs);

    num_atoms = size(a.atoms_data)[1];
    num_clouds = size(a.cloud_data)[1];
    clouds_per_atom = ceil(Int,num_clouds/num_atoms);

    ρa = a.atoms_data[:,4];
    ρb = b.atoms_data[:,4];
    for i in 1:clouds_per_atom
        new_ρ = a.cloud_data[i:clouds_per_atom:end,4];
        new_ρ .*= a.cloud_data[i:clouds_per_atom:end,6];
        ρa -= new_ρ;

        new_ρ = b.cloud_data[i:clouds_per_atom:end,4];
        new_ρ .*= b.cloud_data[i:clouds_per_atom:end,6];
        ρb -= new_ρ;
    end

    return vcat(ρa,ρb);
end

function ChargesAtDistanceOO_ECP(dist::Float64)
    # Returns the energy of the uncorelated model of two molecules of water
    # on the OH reaction coordinate at distance dist.
    return ChargesAtDistanceOO(dist,model1_xc_coeffs_ECP_Pol,"H2O_ecp_fitted_data.txt");
end

function ChargesAtDistanceOO_FullE(dist::Float64)
    # Returns the energy of the uncorelated model of two molecules of water
    # on the OH reaction coordinate at distance dist.
    return ChargesAtDistanceOO(dist,model1_xc_coeffs_FullE_Pol,"H2O_FullE_fitted_data.txt");
end

function ChargesAtDistanceOH(dist::Float64,pol_xc_coeffs::Vector{Float64},
    density_file::String)
    # Returns the energy of the uncorelated model of two molecules of water
    # on the OH reaction coordinate at distance dist.
    a = ReadMolecule(density_file,Float64);
    b = copy(a);
    
    dr = (a.atoms_data[2,1:3] - a.atoms_data[1,1:3]);
    OH_dist = sqrt(sum(dr.^2.0));
    dr += dr.*(dist/OH_dist);

    disp_and_angs = zeros(Float64,6);
    disp_and_angs[4:6] = dr;
    disp_and_angs[3] = bond_θ;

    MoveAndRotateMolec!(b,disp_and_angs);
    PolarizeMolecules!([a,b],pol_xc_coeffs);

    num_atoms = size(a.atoms_data)[1];
    num_clouds = size(a.cloud_data)[1];
    clouds_per_atom = ceil(Int,num_clouds/num_atoms);

    r2 = a.atoms_data[2,1:3];
    r3 = a.atoms_data[3,1:3];
    r4 = b.atoms_data[1,1:3];

    ρa = a.atoms_data[:,4];
    ρb = b.atoms_data[:,4];
    for i in 1:clouds_per_atom
        new_ρ = a.cloud_data[i:clouds_per_atom:end,4];
        new_ρ .*= a.cloud_data[i:clouds_per_atom:end,6];
        ρa -= new_ρ;

        new_ρ = b.cloud_data[i:clouds_per_atom:end,4];
        new_ρ .*= b.cloud_data[i:clouds_per_atom:end,6];
        ρb -= new_ρ;
    end

    return vcat(ρa,ρb);
end

function ChargesAtDistanceOH_ECP(dist::Float64)
    # Returns the energy of the uncorelated model of two molecules of water
    # on the OH reaction coordinate at distance dist.
    # density_file = "H2O_ecp_fitted_data.txt";
    return ChargesAtDistanceOH(dist,model1_xc_coeffs_ECP_Pol,
        "H2O_ecp_fitted_data.txt");
end

function ChargesAtDistanceOH_FullE(dist::Float64)
    # Returns the energy of the uncorelated model of two molecules of water
    # on the OH reaction coordinate at distance dist.
    # density_file = "H2O_ecp_fitted_data.txt";
    return ChargesAtDistanceOH(dist,model1_xc_coeffs_FullE_Pol,
        "H2O_FullE_fitted_data.txt");
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

react_coord = zeros(Float64,0);
ecp_mull_charges = zeros(Float64,0,6);
fulle_mull_charges = zeros(Float64,0,6);

model1_xc_coeffs_ECP = ReadXCCoeffs("ECP_XCCoeffs.txt");
model1_xc_coeffs_FullE = ReadXCCoeffs("FullE_XCCoeffs.txt");
model1_xc_coeffs_ECP_Pol = ReadXCCoeffs("ECP_XCCoeffs_Pol.txt");
model1_xc_coeffs_FullE_Pol = ReadXCCoeffs("FullE_XCCoeffs_Pol.txt");

ECP_order = ceil(Int,length(model1_xc_coeffs_ECP)/4.0) - 1;
FullE_order = ceil(Int,length(model1_xc_coeffs_FullE)/4.0) - 1;

bond_θ = -(π/2.0)*(13.0/25.0);
str_react_coord = "OO";

for i in 0:199
    if str_react_coord == "OH"
        base_name = "Training Data/Gaussian Data/FullE/OH Coord/H2O_H2O_";
        local file_name = base_name * string(13+i*25) * ".log";
        new_charges = ReadGaussianMullikenCharges(file_name)';
        global fulle_mull_charges = vcat(fulle_mull_charges,new_charges);
        mol_geom = GetGeometryFromLogFile(file_name);
        aux_geom = mol_geom[2,:] - mol_geom[4,:];
        aux_react_coord = sqrt(sum(aux_geom.^2.0));
        append!(react_coord,aux_react_coord);
    else
        base_name = "Training Data/Gaussian Data/FullE/OO Coord/H2O_H2O_";
        local file_name = base_name * string(i*25) * ".log" 
        mol_geom = GetGeometryFromLogFile(file_name);
        new_charges = ReadGaussianMullikenCharges(file_name)';
        global fulle_mull_charges = vcat(fulle_mull_charges,new_charges);
        aux_geom = mol_geom[1,:] - mol_geom[4,:];
        aux_react_coord = sqrt(sum(aux_geom.^2.0));
        append!(react_coord,aux_react_coord);
    end
end

for i in 0:199
    if str_react_coord == "OH"
        base_name = "Training Data/Gaussian Data/ECP/OH Coord/H2O_H2O_";
        local file_name = base_name * string(13+i*25) * ".log";
        new_charges = ReadGaussianMullikenCharges(file_name)';
        global ecp_mull_charges = vcat(ecp_mull_charges,new_charges);
        mol_geom = GetGeometryFromLogFile(file_name);
        aux_geom = mol_geom[2,:] - mol_geom[4,:];
    else
        base_name = "Training Data/Gaussian Data/ECP/OO Coord/H2O_H2O_";
        local file_name = base_name * string(i*25) * ".log" 
        mol_geom = GetGeometryFromLogFile(file_name);
        new_charges = ReadGaussianMullikenCharges(file_name)';
        global ecp_mull_charges = vcat(ecp_mull_charges,new_charges);
        aux_geom = mol_geom[1,:] - mol_geom[4,:];
    end
end

react_coord ./=  0.52917721090380;

k = 15;
ecp_mull_charge = ecp_mull_charges[k,:];
fulle_mull_charge = fulle_mull_charges[k,:];

if str_react_coord == "OH"
    model_ecp_charge = ChargesAtDistanceOH_ECP(react_coord[k]);
    model_fulle_charge = ChargesAtDistanceOH_FullE(react_coord[k]);
else
    model_ecp_charge = ChargesAtDistanceOO_ECP(react_coord[k]);
    model_fulle_charge = ChargesAtDistanceOO_FullE(react_coord[k]);
end

p1 = groupedbar([fulle_mull_charge model_fulle_charge ecp_mull_charge model_ecp_charge],
label=["CCSD(T)/cc-pVTZ" "This work (Full E. fit)" "CCSD(T)/CEP-31G" "This work (ECP fit)"]);
plot!(legend=:outertop,legendcolumns=2);
plot!(ylims=(-1,0.5));
plot!(xticks=(1:6,[]));
plot!(ylabel=L"$\rho$");
aux_text = L"$\Delta L = "*string(round(react_coord[k]*1000)/1000)*L"\  \AA$"; 
annotate!(6.5,-0.65,text(aux_text,:center,:right,8));

k = 30;
ecp_mull_charge = ecp_mull_charges[k,:];
fulle_mull_charge = fulle_mull_charges[k,:];

if str_react_coord == "OH"
    model_ecp_charge = ChargesAtDistanceOH_ECP(react_coord[k]);
    model_fulle_charge = ChargesAtDistanceOH_FullE(react_coord[k]);
else
    model_ecp_charge = ChargesAtDistanceOO_ECP(react_coord[k]);
    model_fulle_charge = ChargesAtDistanceOO_FullE(react_coord[k]);
end

p2 = groupedbar([fulle_mull_charge model_fulle_charge ecp_mull_charge model_ecp_charge]);
plot!(legend=false);
plot!(ylims=(-1,0.5));
plot!(xticks=(1:6,[]));
plot!(ylabel=L"$\rho$");
aux_text = L"$\Delta L = "*string(round(react_coord[k]*1000)/1000)*L"\  \AA$"; 
annotate!(6.5,-0.65,text(aux_text,:center,:right,8));

k = 60;
ecp_mull_charge = ecp_mull_charges[k,:];
fulle_mull_charge = fulle_mull_charges[k,:];

if str_react_coord == "OH"
    model_ecp_charge = ChargesAtDistanceOH_ECP(react_coord[k]);
    model_fulle_charge = ChargesAtDistanceOH_FullE(react_coord[k]);
else
    model_ecp_charge = ChargesAtDistanceOO_ECP(react_coord[k]);
    model_fulle_charge = ChargesAtDistanceOO_FullE(react_coord[k]);
end

p3 = groupedbar([fulle_mull_charge model_fulle_charge ecp_mull_charge model_ecp_charge]);
plot!(legend=false);
plot!(ylims=(-1,0.5));
plot!(xlabel="Atom Number");
plot!(xguidefontsize=10);
plot!(ylabel=L"$\rho$");
aux_text = L"$\Delta L = "*string(round(react_coord[k]*1000)/1000)*L"\  \AA$"; 
annotate!(6.5,-0.65,text(aux_text,:center,:right,8));

s1 = 0.45;
s2 = (1-s1)*(0.55);
s3 = (1-s1)*(0.45);

plot(p1,p2,p3,layout=grid(3,1,heights=(s1,s2,s3)))

PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.DEFAULT_PREAMBLE; "\\usepackage{amsmath,accents}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage{accents}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage[utf8]{inputenc}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage[T1]{fontenc}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage{fontspec}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage{mathptmx}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage[fontsize=10pt]{fontsize}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\pgfplotsset{legend image code/.code={ \\draw [#1] (0cm,-0.1cm) rectangle (0.6cm,0.1cm);},}"];

# PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.DEFAULT_PREAMBLE; "\\usepackage[bitstream-charter]{mathdesign}"];
# PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage{amsmath,accents}"];
# PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage{accents}"];
# PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage[utf8]{inputenc}"];
# PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage[T1]{fontenc}"];
# PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage{fontspec}"];
# PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage[fontsize=10pt]{fontsize}"];
# PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\pgfplotsset{legend image code/.code={ \\draw [#1] (0cm,-0.1cm) rectangle (0.6cm,0.1cm);},}"];

model_OH_ecp_charges_15 = ChargesAtDistanceOH_ECP(react_coord[15]);
model_OH_fulle_charges_15 = ChargesAtDistanceOH_FullE(react_coord[15]);
model_OO_ecp_charges_5 = ChargesAtDistanceOO_ECP(react_coord[5]);
model_OO_fulle_charges_5 = ChargesAtDistanceOO_FullE(react_coord[5]);

model_OH_ecp_charges_30 = ChargesAtDistanceOH_ECP(react_coord[30]);
model_OH_fulle_charges_30 = ChargesAtDistanceOH_FullE(react_coord[30]);
model_OO_ecp_charges_25 = ChargesAtDistanceOO_ECP(react_coord[25]);
model_OO_fulle_charges_25 = ChargesAtDistanceOO_FullE(react_coord[25]);

model_OH_ecp_charges_60 = ChargesAtDistanceOH_ECP(react_coord[60]);
model_OH_fulle_charges_60 = ChargesAtDistanceOH_FullE(react_coord[60]);
model_OO_ecp_charges_55 = ChargesAtDistanceOO_ECP(react_coord[55]);
model_OO_fulle_charges_55 = ChargesAtDistanceOO_FullE(react_coord[55]);

if str_react_coord == "OH"
    F1 = @pgf Axis(
        {
            "axis line style={draw=none}",
            "tick style={draw=none}",
            "xticklabels = \\empty",
            "yticklabels = \\empty",
            "draw = none",
            clip = false,
            xmin = 0.0,
            xmax = 1.0,
            ymin = 0.0,
            ymax = 1.0,
            width = 210, 
            height = 45,
        },
        Plot(
            {
                no_marks,
            },
            Coordinates([0],[0]),
            "node[] at (0.5,0.5) {\\includegraphics{/Users/joseantoniorodriguesromero/Documents/GitHub/PolElectronDensityForceModel/labeled_H2O_OH_diag.pdf}}"
        ),
    )
    fp1 = @pgf Axis(
        {
            ybar = 0,
            "bar width = 5.5pt",
            legend_style =
            {
                at = Coordinate(1.0 + (30/230)/2.0, 2.5),
                anchor = "north",
                anchor = "center",
                legend_columns = 2
            },
            ylabel = raw"$\left \lvert \Delta q \right \rvert$",
            ymin = 0.0,
            ymax = 0.2,
            xmin = 0.25,
            xmax = 6.75,
            width = 230, 
            height = 85,
            xtick = "1,2,3,4,5,6",
            xticklabels = "\\empty",
            ytick = "0.0,0.1,0.2",
            yticklabels = L"$0.0$,$0.1$,$0.2$",
            "xtick align = inside",
            "grid = both",
            "grid style={line width=.1pt, draw=gray!10},major grid style={line width=.2pt,draw=gray!25}",
        },
        Plot(
            {fill = theme_palette(:auto)[1]},
            Coordinates(
                [
                    (1, abs(fulle_mull_charges[15,1] - model_OH_fulle_charges_15[1])), 
                    (2, abs(fulle_mull_charges[15,2] - model_OH_fulle_charges_15[2])), 
                    (3, abs(fulle_mull_charges[15,3] - model_OH_fulle_charges_15[3])), 
                    (4, abs(fulle_mull_charges[15,4] - model_OH_fulle_charges_15[4])), 
                    (5, abs(fulle_mull_charges[15,5] - model_OH_fulle_charges_15[5])), 
                    (6, abs(fulle_mull_charges[15,6] - model_OH_fulle_charges_15[6]))  
                ]
            )
        ),
        Plot(
            {fill = theme_palette(:auto)[2]},
            Coordinates(
                [
                    (1, abs(ecp_mull_charges[15,1] - model_OH_ecp_charges_15[1])), 
                    (2, abs(ecp_mull_charges[15,2] - model_OH_ecp_charges_15[2])), 
                    (3, abs(ecp_mull_charges[15,3] - model_OH_ecp_charges_15[3])), 
                    (4, abs(ecp_mull_charges[15,4] - model_OH_ecp_charges_15[4])), 
                    (5, abs(ecp_mull_charges[15,5] - model_OH_ecp_charges_15[5])), 
                    (6, abs(ecp_mull_charges[15,6] - model_OH_ecp_charges_15[6]))  
                ]
            ),
            "node[] at (5.75,0.15) {\\fontsize{8pt}{8pt}\\selectfont \$\\Delta L = "*string(round(react_coord[15],digits=3))*"\\ \\mathring{\\mathrm{A}}\$}",
        ),
        LegendEntry({anchor = "west"},L"\fontsize{8pt}{8pt}\selectfont $\mathrm{Full \ electron \ variant \hspace{0.5cm}}$"),
        LegendEntry({anchor = "west"},L"\fontsize{8pt}{8pt}\selectfont $\mathrm{ECP \ variant}$"),
        "\\pgfplotsset{legend image code/.code={ \\draw [#1] (0cm,-0.1cm) rectangle (0.6cm,0.1cm);},}",
    )

    fp2 = @pgf Axis(
        {
            ybar = 0,
            "bar width = 5.5pt",
            ylabel = raw"$\left \lvert \Delta q \right \rvert$",
            ymin = 0.0,
            ymax = 0.2,
            xmin = 0.25,
            xmax = 6.75,
            width = 230, 
            height = 85,
            xtick = "1,2,3,4,5,6",
            xticklabels = "\\empty",
            ytick = "0.0,0.1,0.2",
            yticklabels = L"$0.0$,$0.1$,$0.2$",
            "xtick align = inside",
            "grid = both",
            "grid style={line width=.1pt, draw=gray!10},major grid style={line width=.2pt,draw=gray!25}",
        },
        Plot(
            {fill = theme_palette(:auto)[1]},
            Coordinates(
                [
                    (1, abs(fulle_mull_charges[30,1] - model_OH_fulle_charges_30[1])), 
                    (2, abs(fulle_mull_charges[30,2] - model_OH_fulle_charges_30[2])), 
                    (3, abs(fulle_mull_charges[30,3] - model_OH_fulle_charges_30[3])), 
                    (4, abs(fulle_mull_charges[30,4] - model_OH_fulle_charges_30[4])), 
                    (5, abs(fulle_mull_charges[30,5] - model_OH_fulle_charges_30[5])), 
                    (6, abs(fulle_mull_charges[30,6] - model_OH_fulle_charges_30[6]))  
                ]
            )
        ),
        Plot(
            {fill = theme_palette(:auto)[2]},
            Coordinates(
                [
                    (1, abs(ecp_mull_charges[30,1] - model_OH_ecp_charges_30[1])), 
                    (2, abs(ecp_mull_charges[30,2] - model_OH_ecp_charges_30[2])), 
                    (3, abs(ecp_mull_charges[30,3] - model_OH_ecp_charges_30[3])), 
                    (4, abs(ecp_mull_charges[30,4] - model_OH_ecp_charges_30[4])), 
                    (5, abs(ecp_mull_charges[30,5] - model_OH_ecp_charges_30[5])), 
                    (6, abs(ecp_mull_charges[30,6] - model_OH_ecp_charges_30[6]))  
                ]
            ),
            "node[] at (5.75,0.15) {\\fontsize{8pt}{8pt}\\selectfont \$\\Delta L = "*string(round(react_coord[30],digits=3))*"\\ \\mathring{\\mathrm{A}}\$}",
        ),
    )

    fp3 = @pgf Axis(
        {
            ybar = 0,
            "bar width = 5.5pt",
            ylabel = raw"$\left \lvert \Delta q \right \rvert$",
            xlabel = L"$\mathrm{Atom \ Number}$",
            ymin = 0.0,
            ymax = 0.2,
            xmin = 0.25,
            xmax = 6.75,
            width = 230, 
            height = 85,
            xtick = "1,2,3,4,5,6",
            xticklabels =  L"$1$,$2$,$3$,$4$,$5$,$6$",
            ytick = "0.0,0.1,0.2",
            yticklabels = L"$0.0$,$0.1$,$0.2$",
            "xtick align = inside",
            "grid = both",
            "grid style={line width=.1pt, draw=gray!10},major grid style={line width=.2pt,draw=gray!25}",
        },
        Plot(
            {fill = theme_palette(:auto)[1]},
            Coordinates(
                [
                    (1, abs(fulle_mull_charges[60,1] - model_OH_fulle_charges_60[1])), 
                    (2, abs(fulle_mull_charges[60,2] - model_OH_fulle_charges_60[2])), 
                    (3, abs(fulle_mull_charges[60,3] - model_OH_fulle_charges_60[3])), 
                    (4, abs(fulle_mull_charges[60,4] - model_OH_fulle_charges_60[4])), 
                    (5, abs(fulle_mull_charges[60,5] - model_OH_fulle_charges_60[5])), 
                    (6, abs(fulle_mull_charges[60,6] - model_OH_fulle_charges_60[6]))  
                ]
            )
        ),
        Plot(
            {fill = theme_palette(:auto)[2]},
            Coordinates(
                [
                    (1, abs(ecp_mull_charges[60,1] - model_OH_ecp_charges_60[1])), 
                    (2, abs(ecp_mull_charges[60,2] - model_OH_ecp_charges_60[2])), 
                    (3, abs(ecp_mull_charges[60,3] - model_OH_ecp_charges_60[3])), 
                    (4, abs(ecp_mull_charges[60,4] - model_OH_ecp_charges_60[4])), 
                    (5, abs(ecp_mull_charges[60,5] - model_OH_ecp_charges_60[5])), 
                    (6, abs(ecp_mull_charges[60,6] - model_OH_ecp_charges_60[6]))  
                ]
            ),
            "node[] at (5.75,0.15) {\\fontsize{8pt}{8pt}\\selectfont \$\\Delta L = "*string(round(react_coord[60],digits=3))*"\\ \\mathring{\\mathrm{A}}\$}",
        ),
    )
else
    F2 = @pgf Axis(
        {
            "axis line style={draw=none}",
            "tick style={draw=none}",
            "xticklabels = \\empty",
            "yticklabels = \\empty",
            "draw = none",
            clip = false,
            xmin = 0.0,
            xmax = 1.0,
            ymin = 0.0,
            ymax = 1.0,
            width = 210, 
            height = 45,
        },
        Plot(
            {
                no_marks,
            },
            Coordinates([0],[0]),
            "node[] at (0.5,0.5) {\\includegraphics{/Users/joseantoniorodriguesromero/Documents/GitHub/PolElectronDensityForceModel/labeled_H2O_OO_diag.pdf}}"
        ),
    )
    fg1 = @pgf Axis(
        {
            ybar = 0,
            "bar width = 5.5pt",
            legend_style =
            {
                at = Coordinate(0.5, 3.0),
                anchor = "north",
                legend_columns = 2
            },
            ymin = 0.0,
            ymax = 0.2,
            xmin = 0.25,
            xmax = 6.75,
            width = 230, 
            height = 85,
            xtick = "1,2,3,4,5,6",
            xticklabels = "\\empty",
            ytick = "0.0,0.1,0.2",
            yticklabels = "\\empty",
            "xtick align = inside",
            "grid = both",
            "grid style={line width=.1pt, draw=gray!10},major grid style={line width=.2pt,draw=gray!25}",
        },
        Plot(
            {fill = theme_palette(:auto)[1]},
            Coordinates(
                [
                    (1, abs(fulle_mull_charges[5,1] - model_OO_fulle_charges_5[1])), 
                    (2, abs(fulle_mull_charges[5,2] - model_OO_fulle_charges_5[2])), 
                    (3, abs(fulle_mull_charges[5,3] - model_OO_fulle_charges_5[3])), 
                    (4, abs(fulle_mull_charges[5,4] - model_OO_fulle_charges_5[4])), 
                    (5, abs(fulle_mull_charges[5,5] - model_OO_fulle_charges_5[5])), 
                    (6, abs(fulle_mull_charges[5,6] - model_OO_fulle_charges_5[6]))  
                ]
            )
        ),
        Plot(
            {fill = theme_palette(:auto)[2]},
            Coordinates(
                [
                    (1, abs(ecp_mull_charges[5,1] - model_OO_ecp_charges_5[1])), 
                    (2, abs(ecp_mull_charges[5,2] - model_OO_ecp_charges_5[2])), 
                    (3, abs(ecp_mull_charges[5,3] - model_OO_ecp_charges_5[3])), 
                    (4, abs(ecp_mull_charges[5,4] - model_OO_ecp_charges_5[4])), 
                    (5, abs(ecp_mull_charges[5,5] - model_OO_ecp_charges_5[5])), 
                    (6, abs(ecp_mull_charges[5,6] - model_OO_ecp_charges_5[6]))  
                ]
            ),
            "node[] at (5.75,0.15) {\\fontsize{8pt}{8pt}\\selectfont \$\\Delta L = "*string(round(react_coord[5],digits=3))*"\\ \\mathring{\\mathrm{A}}\$}",
        ),
        # LegendEntry({anchor = "west"},L"\fontsize{8pt}{8pt}\selectfont $\mathrm{CCSD(T)/cc\text{-}pVTZ}$"),
        # LegendEntry({anchor = "west"},L"\fontsize{8pt}{8pt}\selectfont $\mathrm{This \ Work \ (Full \ E. \ fit)}$"),
        # LegendEntry({anchor = "west"},L"\fontsize{8pt}{8pt}\selectfont $\mathrm{CCSD(T)/CEP\text{-}31G}$ \qquad"),
        # LegendEntry({anchor = "west"},L"\fontsize{8pt}{8pt}\selectfont $\mathrm{This \ Work \ (ECP \ fit)}$"),
        # "\\pgfplotsset{legend image code/.code={ \\draw [#1] (0cm,-0.1cm) rectangle (0.6cm,0.1cm);},}",
    )

    fg2 = @pgf Axis(
        {
            ybar = 0,
            "bar width = 5.5pt",
            ymin = 0.0,
            ymax = 0.2,
            xmin = 0.25,
            xmax = 6.75,
            width = 230, 
            height = 85,
            xtick = "1,2,3,4,5,6",
            xticklabels = "\\empty",
            ytick = "0.0,0.1,0.2",
            yticklabels = "\\empty",
            "xtick align = inside",
            "grid = both",
            "grid style={line width=.1pt, draw=gray!10},major grid style={line width=.2pt,draw=gray!25}",
        },
        Plot(
            {fill = theme_palette(:auto)[1]},
            Coordinates(
                [
                    (1, abs(fulle_mull_charges[25,1] - model_OO_fulle_charges_25[1])), 
                    (2, abs(fulle_mull_charges[25,2] - model_OO_fulle_charges_25[2])), 
                    (3, abs(fulle_mull_charges[25,3] - model_OO_fulle_charges_25[3])), 
                    (4, abs(fulle_mull_charges[25,4] - model_OO_fulle_charges_25[4])), 
                    (5, abs(fulle_mull_charges[25,5] - model_OO_fulle_charges_25[5])), 
                    (6, abs(fulle_mull_charges[25,6] - model_OO_fulle_charges_25[6]))  
                ]
            )
        ),
        Plot(
            {fill = theme_palette(:auto)[2]},
            Coordinates(
                [
                    (1, abs(ecp_mull_charges[25,1] - model_OO_ecp_charges_25[1])), 
                    (2, abs(ecp_mull_charges[25,2] - model_OO_ecp_charges_25[2])), 
                    (3, abs(ecp_mull_charges[25,3] - model_OO_ecp_charges_25[3])), 
                    (4, abs(ecp_mull_charges[25,4] - model_OO_ecp_charges_25[4])), 
                    (5, abs(ecp_mull_charges[25,5] - model_OO_ecp_charges_25[5])), 
                    (6, abs(ecp_mull_charges[25,6] - model_OO_ecp_charges_25[6]))  
                ]
            ),
            "node[] at (5.75,0.15) {\\fontsize{8pt}{8pt}\\selectfont \$\\Delta L = "*string(round(react_coord[25],digits=3))*"\\ \\mathring{\\mathrm{A}}\$}",
        )
    )

    fg3 = @pgf Axis(
        {
            ybar = 0,
            "bar width = 5.5pt",
            xlabel = L"$\mathrm{Atom \ Number}$",
            ymin = 0.0,
            ymax = 0.2,
            xmin = 0.25,
            xmax = 6.75,
            width = 230, 
            height = 85,
            xtick = "1,2,3,4,5,6",
            xticklabels =  L"$1$,$2$,$3$,$4$,$5$,$6$",
            ytick = "0.0,0.1,0.2",
            yticklabels = "\\empty",
            "xtick align = inside",
            "grid = both",
            "grid style={line width=.1pt, draw=gray!10},major grid style={line width=.2pt,draw=gray!25}",
        },
        Plot(
            {fill = theme_palette(:auto)[1]},
            Coordinates(
                [
                    (1, abs(fulle_mull_charges[55,1] - model_OO_fulle_charges_55[1])), 
                    (2, abs(fulle_mull_charges[55,2] - model_OO_fulle_charges_55[2])), 
                    (3, abs(fulle_mull_charges[55,3] - model_OO_fulle_charges_55[3])), 
                    (4, abs(fulle_mull_charges[55,4] - model_OO_fulle_charges_55[4])), 
                    (5, abs(fulle_mull_charges[55,5] - model_OO_fulle_charges_55[5])), 
                    (6, abs(fulle_mull_charges[55,6] - model_OO_fulle_charges_55[6]))  
                ]
            )
        ),
        Plot(
            {fill = theme_palette(:auto)[2]},
            Coordinates(
                [
                    (1, abs(ecp_mull_charges[55,1] - model_OO_ecp_charges_55[1])), 
                    (2, abs(ecp_mull_charges[55,2] - model_OO_ecp_charges_55[2])), 
                    (3, abs(ecp_mull_charges[55,3] - model_OO_ecp_charges_55[3])), 
                    (4, abs(ecp_mull_charges[55,4] - model_OO_ecp_charges_55[4])), 
                    (5, abs(ecp_mull_charges[55,5] - model_OO_ecp_charges_55[5])), 
                    (6, abs(ecp_mull_charges[55,6] - model_OO_ecp_charges_55[6]))  
                ]
            ),
            "node[] at (5.75,0.15) {\\fontsize{8pt}{8pt}\\selectfont \$\\Delta L = "*string(round(react_coord[55],digits=3))*"\\ \\mathring{\\mathrm{A}}\$}",
        ),
    )
end


@pgf GroupPlot({group_style = { group_size = "2 by 4","horizontal sep = 30"}},F1,F2,fp1,fg1,fp2,fg2,fp3,fg3)