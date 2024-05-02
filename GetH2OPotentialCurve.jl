using Plots, SpecialFunctions, LaTeXStrings, PGFPlotsX
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

function EnergyAtDistanceOO(dist::Float64,order::Int,
    pol_xc_coeffs::Vector{Float64},density_file::String)
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

    energies = zeros(Float64,1+4*(order+1));
    energies[1] = NaiveEnergyFromDensity([a,b]);
    energies[2:end] = XCEnergyFromDensity([a,b],order);
    return energies;
end

function EnergyAtDistanceOO_ECP(dist::Float64,order::Int)
    # Returns the energy of the uncorelated model of two molecules of water
    # on the OH reaction coordinate at distance dist.
    return EnergyAtDistanceOO(dist,order,model1_xc_coeffs_ECP_Pol,"H2O_ecp_fitted_data.txt");
end

function EnergyAtDistanceOO_FullE(dist::Float64,order::Int)
    # Returns the energy of the uncorelated model of two molecules of water
    # on the OH reaction coordinate at distance dist.
    return EnergyAtDistanceOO(dist,order,model1_xc_coeffs_FullE_Pol,"H2O_FullE_fitted_data.txt");
end

function EnergyAtDistanceOH(dist::Float64,order::Int,
    pol_xc_coeffs::Vector{Float64},density_file::String)
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

    energies = zeros(Float64,1+4*(order+1),1);
    energies[1] = NaiveEnergyFromDensity([a,b]);
    energies[2:end] = XCEnergyFromDensity([a,b],order);
    return energies;
end

function EnergyAtDistanceOH_ECP(dist::Float64,order::Int)
    # Returns the energy of the uncorelated model of two molecules of water
    # on the OH reaction coordinate at distance dist.
    # density_file = "H2O_ecp_fitted_data.txt";
    return EnergyAtDistanceOH(dist,order,model1_xc_coeffs_ECP_Pol,
        "H2O_ecp_fitted_data.txt");
end

function EnergyAtDistanceOH_FullE(dist::Float64,order::Int)
    # Returns the energy of the uncorelated model of two molecules of water
    # on the OH reaction coordinate at distance dist.
    # density_file = "H2O_ecp_fitted_data.txt";
    return EnergyAtDistanceOH(dist,order,model1_xc_coeffs_FullE_Pol,
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

file_name = "Training Data/Gaussian Data/FullE/OH Coord/H2O.log";
e0 = GetEnergyFromLogFile(file_name);
surf_energy = zeros(Float64,0);
react_coord = zeros(Float64,0);

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
        append!(surf_energy,GetEnergyFromLogFile(file_name));
        mol_geom = GetGeometryFromLogFile(file_name);
        aux_geom = mol_geom[2,:] - mol_geom[4,:];
        aux_react_coord = sqrt(sum(aux_geom.^2.0));
        append!(react_coord,aux_react_coord);
    else
        base_name = "Training Data/Gaussian Data/FullE/OO Coord/H2O_H2O_";
        local file_name = base_name * string(i*25) * ".log" 
        mol_geom = GetGeometryFromLogFile(file_name);
        append!(surf_energy,GetEnergyFromLogFile(file_name));
        aux_geom = mol_geom[1,:] - mol_geom[4,:];
        aux_react_coord = sqrt(sum(aux_geom.^2.0));
        append!(react_coord,aux_react_coord);
    end
end

a0 = 0.529177210903;
surf_energy .-= 2.0*e0;

# model1_coord = (react_coord./a0);
model1_coord = zeros(Float64,801,1);
model1_coord[:] = 1.0 .+ collect(0:800).*(9.0/800.0);

if str_react_coord == "OH"
    t = time();
    ECP_all_model1_e = [];
    FullE_all_model1_e = [];
    for aux_coord in model1_coord
        push!(ECP_all_model1_e,EnergyAtDistanceOH_ECP(aux_coord,ECP_order));
        push!(FullE_all_model1_e,EnergyAtDistanceOH_FullE(aux_coord,FullE_order));
    end
    t = time() - t;

    println("Time taken (this model): "*string(t));
else
    t = time();
    ECP_all_model1_e = [];
    FullE_all_model1_e = [];
    for aux_coord in model1_coord
        push!(ECP_all_model1_e,EnergyAtDistanceOO_ECP(aux_coord,ECP_order));
        push!(FullE_all_model1_e,EnergyAtDistanceOO_FullE(aux_coord,FullE_order));
    end
    t = time() - t;
    println("Time taken (this model): "*string(t));
end

t = time();
for i in eachindex(model1_coord)
    density_file = "H2O_ecp_fitted_data.txt";
    local a = ReadMolecule(density_file,Float64);
    local b = ReadMolecule(density_file,Float64);
    dist = model1_coord[i];

    if str_react_coord == "OH"
        dr = (a.atoms_data[2,1:3] - a.atoms_data[1,1:3]);
        OH_dist = sqrt(sum(dr.^2.0));
        dr += dr.*(dist/OH_dist);

        disp_and_angs = zeros(Float64,6);
        disp_and_angs[4:6] = dr;
        disp_and_angs[3] = bond_θ;

        MoveAndRotateMolec!(b,disp_and_angs);
    else
        disp_and_angs = zeros(Float64,6);
        disp_and_angs[3] = π;
        disp_and_angs[5] = dist;

        MoveAndRotateMolec!(a,disp_and_angs);
    end
end

t = time() - t;
println("Baseline: "*string(t));

ECP_all_model1_e = reduce(hcat,ECP_all_model1_e)';
ECP_model1_e = ECP_all_model1_e[:,1];
ECP_xc_model1 = ECP_all_model1_e[:,2:end];

FullE_all_model1_e = reduce(hcat,FullE_all_model1_e)';
FullE_model1_e = FullE_all_model1_e[:,1];
FullE_xc_model1 = FullE_all_model1_e[:,2:end];

ECP_model1_e += ECP_xc_model1*model1_xc_coeffs_ECP;
FullE_model1_e += FullE_xc_model1*model1_xc_coeffs_FullE;

model1_coord .*= a0;
kjmol = 2625.5002;

ECP_H2O = [ReadMolecule("H2O_ecp_fitted_data.txt")];
PolarizeMolecules!(ECP_H2O,model1_xc_coeffs_ECP_Pol);

ECP_e0 = NaiveEnergyFromDensity(ECP_H2O);
ECP_e0 += dot(XCEnergyFromDensity(ECP_H2O,ECP_order),model1_xc_coeffs_ECP);
ECP_model1_e .-= 2.0*ECP_e0;

FullE_H2O = [ReadMolecule("H2O_fullE_fitted_data.txt")];
PolarizeMolecules!(FullE_H2O,model1_xc_coeffs_FullE_Pol);

FullE_e0 = NaiveEnergyFromDensity(FullE_H2O);
FullE_e0 += dot(XCEnergyFromDensity(FullE_H2O,FullE_order),model1_xc_coeffs_FullE);
FullE_model1_e .-= 2.0*FullE_e0;

ECP_model1_e *= kjmol;
FullE_model1_e *= kjmol;
surf_energy *= kjmol;

l_width = 2.5;
legends = ["CCSD(T)/cc-pVTZ";"This Work (ECP fit)";
"This Work (Full E. fit)";"GAFF";"MMFF94S";"TIP3P";"MB-Pol";"EFP"];

p = plot(model1_coord,ECP_model1_e,labels=legends[2],linewidth=l_width);
plot!(model1_coord,FullE_model1_e,labels=legends[3],linewidth=l_width);
plot!(react_coord,surf_energy,labels=legends[1],linewidth=l_width);

if str_react_coord == "OH"
    plot!(ylims=(-50,750),xlims=(0.75,2));
    plot!(yticks=0:250:1000);

    plot!(legend=false);
    plot!(size=(600,255));
    plot!(left_margin=2Plots.mm);
    plot!(right_margin=2Plots.mm);
    plot!(bottom_margin=2Plots.mm,top_margin=2Plots.mm);

    plot!(ylabel=L"\Delta E \ \ [\textrm{kJ/mol}]");
    plot!(xlabel=L"\Delta L \ \ [\AA]");

    p_copy = deepcopy(p);

    plot!(size=(600,400));
    plot!(yticks=-30:15:30);
    plot!(ylims=(-30,30),xlims=(0.5,5));

    # plot!(yticks=([0.015,0,-0.015],[]))
    # plot!(xticks=(1:4,[]));
    # plot!(xlabel="");

    plot!(xlabel="");
    plot!(legend=:outertop,legendcolumns=3);
    # plot!(legend=:outertop, legendguide=:bottomright)
    # plot!(top_margin=5Plots.mm)

    p = plot(p,p_copy,layout=grid(2,1,heights=(4.5/8,3.5/8)))
else
    plot!(ylims=(-50,800),xlims=(1.0,3.5));
    plot!(yticks=0:200:800);

    plot!(size=(600,255));
    plot!(legend=:outertop,legendcolumns=3);
    plot!(left_margin=2Plots.mm);
    plot!(right_margin=2Plots.mm);
    plot!(bottom_margin=2Plots.mm,top_margin=2Plots.mm);

    # plot!(ylims=(-50,600),xlims=(0,5));
    # plot!(yticks=0:150:600);

    # plot!(yticks=([0.015,0,-0.015],[]))

    plot!(xlabel=L"\Delta L \ \ [\AA]");
    plot!(ylabel=L"\Delta E \ \ [\textrm{kJ/mol}]");

    # plot!(xticks=(1:4,[]));
    # plot!(xlabel="");

    plot!(size=(600,255));
    plot!(legend=:outertop,legendcolumns=3);
    plot!(left_margin=2Plots.mm);
    plot!(right_margin=2Plots.mm);
    plot!(bottom_margin=2Plots.mm,top_margin=2Plots.mm)
    # plot!(legend=:outertop, legendguide=:bottomright)
    # plot!(top_margin=5Plots.mm)
end

PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.DEFAULT_PREAMBLE; "\\usepackage{amsmath,accents}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage{accents}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage[utf8]{inputenc}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage[T1]{fontenc}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage{fontspec}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage{mathptmx}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage[fontsize=10pt]{fontsize}"];

# PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.DEFAULT_PREAMBLE; "\\usepackage[bitstream-charter]{mathdesign}"];
# PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage{amsmath,accents}"];
# PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage{accents}"];
# PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage[utf8]{inputenc}"];
# PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage[T1]{fontenc}"];
# PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage{fontspec}"];
# PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage[fontsize=10pt]{fontsize}"];

aux_cond_ecp = (abs.(ECP_model1_e) .< 4000);
aux_cond_fulle = (abs.(FullE_model1_e) .< 4000);

if str_react_coord == "OO"
    fp1 = @pgf Axis(
        {
            xmajorgrids,
            ymajorgrids,
            legend_style =
            {
                at = Coordinate(0.5,-0.5 - 1.5*(30.0/125.0)),
                anchor = "center",
                legend_columns = 1
            },
            xmin = 1.0,
            xmax = 3.5,
            ymin = -50,
            ymax = 800,
            xtick = "1.0,1.5,2.0,2.5,3.0,3.5",
            ytick = "0,200,400,600,800",
            xticklabels = L"$1.0$,$1.5$,$2.0$,$2.5$,$3.0$,$3.5$",
            xlabel = L"$\Delta L \ [ \mathring{\mathrm{A}} ]$",
            # ylabel = L"$\Delta E \ \left[ \mathrm{kJ/mol} \right]$",
            width = 210, 
            height = 125,
            "grid style={line width=.1pt, draw=gray!10},major grid style={line width=.2pt,draw=gray!25}",
        },
        Plot(
            {
                color = theme_palette(:auto)[1],
                style = {"thick"},
                no_marks,
            },
            Coordinates(model1_coord[aux_cond_ecp][:],
                ECP_model1_e[aux_cond_ecp][:]),
        ),
        Plot(
            {
                color = theme_palette(:auto)[2],
                no_marks,
                style = {"thick"},
            },
            Coordinates(model1_coord[aux_cond_fulle][:],
                FullE_model1_e[aux_cond_fulle][:]),
        ),
        Plot(
            {
                color = theme_palette(:auto)[3],
                no_marks,
                style = {"thick"},
            },
            Coordinates(react_coord,surf_energy),
        ),
        LegendEntry({anchor = "west"},L"\fontsize{8pt}{8pt}\selectfont $\mathrm{This \ Work \ (ECP \ variant)}$"),
        LegendEntry({anchor = "west"},L"\fontsize{8pt}{8pt}\selectfont $\mathrm{This \ Work \ (Full \ electron \ variant)}$"),
        LegendEntry({anchor = "west"},L"\fontsize{8pt}{8pt}\selectfont $\mathrm{CCSD(T)/cc\text{-}pVTZ}$")
    )
elseif str_react_coord == "OH"
    gp1 = @pgf Axis(
        {
            xmajorgrids,
            ymajorgrids,
            xmin = 0.5,
            xmax = 5,
            ymin = -30,
            ymax = 30,
            xlabel = L"$\Delta L \ [ \mathring{\mathrm{A}} ]$",
            xtick = "1,2,3,4,5",
            ytick = "-30,-15,0,15,30",
            ylabel = L"$\Delta E \ \left[ \mathrm{kJ/mol} \right]$",
            width = 210, 
            height = 125,
            "grid style={line width=.1pt, draw=gray!10},major grid style={line width=.2pt,draw=gray!25}",
        },
        Plot(
            {
                color = theme_palette(:auto)[1],
                style = {"thick"},
                no_marks,
            },
            Coordinates(model1_coord[aux_cond_ecp][:],
                ECP_model1_e[aux_cond_ecp][:]),
        ),
        Plot(
            {
                color = theme_palette(:auto)[2],
                no_marks,
                style = {"thick"},
            },
            Coordinates(model1_coord[aux_cond_fulle][:],
                FullE_model1_e[aux_cond_fulle][:]),
        ),
        Plot(
            {
                color = theme_palette(:auto)[3],
                no_marks,
                style = {"thick"},
            },
            Coordinates(react_coord,surf_energy),
        )
    )
    gp2 = @pgf Axis(
        {
            xmajorgrids,
            ymajorgrids,
            xmin = 0.75,
            xmax = 2.0,
            ymin = -50,
            ymax = 750,
            xtick = "0.75,1.0,1.25,1.5,1.75,2.0",
            xticklabels = L"$0.75$,$1.00$,$1.25$,$1.50$,$1.75$,$2.00$",
            ytick = "0,250,500,750",
            ylabel = L"$\Delta E \ \left[ \mathrm{kJ/mol} \right]$",
            width = 210, 
            height = 125,
            "grid style={line width=.1pt, draw=gray!10},major grid style={line width=.2pt,draw=gray!25}",
        },
        Plot(
            {
                color = theme_palette(:auto)[1],
                style = {"thick"},
                no_marks,
            },
            Coordinates(model1_coord[aux_cond_ecp][:],
                ECP_model1_e[aux_cond_ecp][:]),
        ),
        Plot(
            {
                color = theme_palette(:auto)[2],
                no_marks,
                style = {"thick"},
            },
            Coordinates(model1_coord[aux_cond_fulle][:],
                FullE_model1_e[aux_cond_fulle][:]),
        ),
        Plot(
            {
                color = theme_palette(:auto)[3],
                no_marks,
                style = {"thick"},
            },
            Coordinates(react_coord,surf_energy),
        ),
        # LegendEntry({anchor = "west"},L"\fontsize{8pt}{8pt}\selectfont $\mathrm{This \ Work \ (ECP \ fit)}$"),
        # LegendEntry({anchor = "west"},L"\fontsize{8pt}{8pt}\selectfont $\mathrm{This \ Work \ (Full \ E. \ fit)}$"),
        # LegendEntry({anchor = "west"},L"\fontsize{8pt}{8pt}\selectfont $\mathrm{CCSD(T)/cc\text{-}pVTZ}$")
    )

    # @pgf GroupPlot({group_style = { group_size = "1 by 2",}},fp2,fp1);
end

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
        "node[] at (0.5,0.5) {\\includegraphics{/Users/joseantoniorodriguesromero/Documents/GitHub/PolElectronDensityForceModel/H2O_OH_diag.pdf}}"
    ),
)

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
        "node[] at (0.5,0.5) {\\includegraphics{/Users/joseantoniorodriguesromero/Documents/GitHub/PolElectronDensityForceModel/H2O_OO_diag.pdf}}"
    ),
)

@pgf GroupPlot({group_style = { group_size = "2 by 3","vertical sep = 30","horizontal sep = 40"}},F1,F2,gp2,fp1,gp1,)