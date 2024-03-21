using Plots, SpecialFunctions, LaTeXStrings, Printf, PGFPlotsX
include("GetPotentialFromAngles.jl")

kjmol = 2625.5002;
a0 = 0.529177210903;
kCalMol_to_Hartree = 0.0015936011;

# PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.DEFAULT_PREAMBLE; "\\usepackage[bitstream-charter]{mathdesign}"];
# PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage{amsmath,mathtools,accents}"];
# PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage{accents}"];
# PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage[utf8]{inputenc}"];
# PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage[T1]{fontenc}"];
# PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage{fontspec}"];
# PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage[fontsize=10pt]{fontsize}"];

PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.DEFAULT_PREAMBLE; "\\usepackage{amsmath,accents}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage{accents}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage[utf8]{inputenc}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage[T1]{fontenc}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage{fontspec}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage{mathptmx}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage[fontsize=10pt]{fontsize}"];

function ReadXCCoeffs(file_name::String)
    # Reads the XC coefficients from the specified text file. The file is 
    # assumed to be in a single column format.
    ret = zeros(Float64,0);
    for line in readlines(file_name)
        push!(ret,parse(Float64,line));
    end

    return ret;
end

xc_coeffs_ECP = ReadXCCoeffs("ECP_XCCoeffs.txt");
xc_coeffs_FullE = ReadXCCoeffs("FullE_XCCoeffs.txt");
xc_coeffs_ECP_Pol = ReadXCCoeffs("ECP_XCCoeffs_Pol.txt");
xc_coeffs_FullE_Pol = ReadXCCoeffs("FullE_XCCoeffs_Pol.txt");

# GAFF Data
N2_N2_GAFF_data = zeros(Float64,0);
for line in readlines("Test Data/GAFF Data/Neutral/N2+N2/GAFF_energies.txt")
    if length(split(line)) > 3
        push!(N2_N2_GAFF_data,parse(Float64,split(line)[end-1]));
    end
end

dist_N2_N2_GAFF_data = 1.0 .+ ((4.0/100.0).*collect(0:100));
N2_N2_GAFF_data = N2_N2_GAFF_data[2:end] .- 2.0*N2_N2_GAFF_data[1];

# MMFF94S Data
N2_N2_MMFF94S_data = zeros(Float64,0);
for line in readlines("Test Data/MMFF94S Data/Neutral/N2+N2/MMFF94S_energies.txt")
    if length(split(line)) > 3
        push!(N2_N2_MMFF94S_data,parse(Float64,split(line)[end-1]));
    end
end

dist_N2_N2_MMFF94S_data = 1.0 .+ ((4.0/100.0).*collect(0:100));
N2_N2_MMFF94S_data = N2_N2_MMFF94S_data[2:end] .- 2.0*N2_N2_MMFF94S_data[1];
N2_N2_MMFF94S_data .*= kjmol*kCalMol_to_Hartree*3.5;

# Gaussian Data
function ReadGaussianEnergy(file_name::String)
    for line in readlines(file_name)
        if contains(line,"CCSD(T)=")
            return parse(Float64,replace(split(line)[end],"D"=>"E"))
        end
    end

    return parse(Float64,"inf");
end

dist_N2_N2_CCSDT_data = 1.0 .+ ((4.0/100.0).*collect(0:100));
N2_N2_CCSDT_data = zeros(Float64,303);
for i in 0:302
    local file_name = "Test Data/Gaussian Data/Neutral/N2+N2/N2_N2_"*string(i)*".log"
    N2_N2_CCSDT_data[i+1] = ReadGaussianEnergy(file_name)
end

CCSDT_e0 = ReadGaussianEnergy("Test Data/Gaussian Data/Neutral/N2+N2/N2.log");
N2_N2_CCSDT_data .-= 2.0*CCSDT_e0;
N2_N2_CCSDT_data .*= kjmol;

# θ = 0
N2_N2_ecp_model_data = zeros(Float64,303);
dist_ecp_model = (1.0 .+ ((4.0/100.0).*collect(0:100)));
model_ecp_data = zeros(Float64,0);
mol_ecp_a = ReadMolecule("N2_ecp_fitted_data.txt");

PolarizeMolecules!([mol_ecp_a],xc_coeffs_ECP_Pol);
order = ceil(Int,length(xc_coeffs_ECP) / 4) - 1;

aux_energy = NaiveEnergyFromDensity([mol_ecp_a]);
aux_energy += (XCEnergyFromDensity([mol_ecp_a],order)*xc_coeffs_ECP)[1];
model_ecp_e0 = 2.0*aux_energy;

for dist in dist_ecp_model
    local mol_ecp_b = copy(mol_ecp_a);
    mol_len = abs(mol_ecp_b.atoms_data[1,1] - mol_ecp_b.atoms_data[2,1]);
    MoveAndRotateMolec!(mol_ecp_b,[0,0,0,dist/a0 + mol_len,0,0]);
    PolarizeMolecules!([mol_ecp_a,mol_ecp_b],xc_coeffs_ECP_Pol);

    local aux_energy = NaiveEnergyFromDensity([mol_ecp_a,mol_ecp_b]);
    local aux_energy += (XCEnergyFromDensity([mol_ecp_a,mol_ecp_b],order)*xc_coeffs_ECP)[1];
    push!(model_ecp_data,(aux_energy-model_ecp_e0)*kjmol*3.0);
end
N2_N2_ecp_model_data[1:3:end] = model_ecp_data;

p1 = plot(dist_N2_N2_GAFF_data,N2_N2_GAFF_data[1:3:end],label="GAFF");
plot!(dist_N2_N2_MMFF94S_data,N2_N2_MMFF94S_data[1:3:end],label="MMFF94S (×3.0)");
plot!(dist_N2_N2_CCSDT_data,N2_N2_CCSDT_data[1:3:end],label="CCSD(T)/cc-pVTZ");
plot!(dist_ecp_model,model_ecp_data,label="This Work (ECP Fit) (×3.0)");
plot!(ylims=(-4,2));
plot!(xlims=(2.5,5));
plot!(yticks=-4:2:2);
plot!(xticks=(2:5,[]));
plot!(ylabel=L"$\Delta E \ \left[ \mathrm{kJ/mol} \right]$");
plot!(legend=:outertop,legendcolumns=2);
annotate!(4.7,-3,text(L"\theta = 0",:center,:left,10));

fg1 = @pgf Axis(
    {
        legend_style =
        {
            at = Coordinate(0.5,1.5),
            anchor = "center",
            legend_columns = 2
        },
        xmajorgrids,
        ymajorgrids,
        width = 280, 
        height = 100,
        xmin = 2.75,
        xmax = 5,
        ymin = -4,
        ymax = 2,
        xtick = "3,3.5,4,4.5,5",
        xticklabels = "\\empty",
        ylabel = L"\Delta E \ \left[ \mathrm{kJ/mol} \right]",
        "grid style={line width=.1pt, draw=gray!10},major grid style={line width=.2pt,draw=gray!25}",
    },
    Plot(
        {
            "no marks",
            style = {"thick"},
            color = theme_palette(:auto)[1],
        },
        Coordinates(dist_ecp_model[model_ecp_data .< 10],model_ecp_data[model_ecp_data .< 10]),
        "node[anchor=center] at (4.75,-3) {\\fontsize{8pt}{8pt}\\selectfont \\color{black} \$ \\theta = 0 \$}",
    ),
    Plot(
        {
            "no marks",
            style = {"thick"},
            color = theme_palette(:auto)[3],
        },
        Coordinates(dist_N2_N2_MMFF94S_data[N2_N2_MMFF94S_data[1:3:end] .< 5.0],N2_N2_MMFF94S_data[1:3:end][N2_N2_MMFF94S_data[1:3:end] .< 5.0]),
    ),
    Plot(
        {
            "no marks",
            style = {"thick"},
            color = theme_palette(:auto)[2],
        },
        Coordinates(dist_N2_N2_CCSDT_data,N2_N2_CCSDT_data[1:3:end]),
    ),
    Plot(
        {
            "no marks",
            style = {"thick"},
            color = theme_palette(:auto)[4],
        },
        Coordinates(dist_N2_N2_GAFF_data[N2_N2_GAFF_data[1:3:end] .< 5.0],N2_N2_GAFF_data[1:3:end][N2_N2_GAFF_data[1:3:end] .< 5.0]),
    ),
    LegendEntry({anchor = "west"},L"\fontsize{8pt}{8pt}\selectfont $\mathrm{This \ Work \ (ECP \ fit)} \ \left( \times 3.0 \right)$"),
    LegendEntry({anchor = "west"},L"\fontsize{8pt}{8pt}\selectfont $\mathrm{MMFF94S} \ \left( \times 3.0 \right)$"),
    LegendEntry({anchor = "west"},L"\fontsize{8pt}{8pt}\selectfont $\mathrm{CCSD(T)/cc\text{-}pVTZ}$"),
    LegendEntry({anchor = "west"},L"\fontsize{8pt}{8pt}\selectfont $\mathrm{GAFF}$"),
);

# θ = π/4
model_ecp_data = zeros(Float64,0);
for dist in dist_ecp_model
    local mol_ecp_b = copy(mol_ecp_a);
    mol_len = abs(mol_ecp_b.atoms_data[1,1] - mol_ecp_b.atoms_data[2,1]);
    MoveAndRotateMolec!(mol_ecp_b,[0,0,0,dist/a0 + mol_len,0,0]);
    MoveAndRotateMolec!(mol_ecp_b,[0,0,π/4,0,0,0]);
    PolarizeMolecules!([mol_ecp_a,mol_ecp_b],xc_coeffs_ECP_Pol);

    local aux_energy = NaiveEnergyFromDensity([mol_ecp_a,mol_ecp_b]);
    local aux_energy += (XCEnergyFromDensity([mol_ecp_a,mol_ecp_b],order)*xc_coeffs_ECP)[1];
    push!(model_ecp_data,(aux_energy-model_ecp_e0)*kjmol*3.0);
end
N2_N2_ecp_model_data[2:3:end] = model_ecp_data;

p2 = plot(dist_N2_N2_GAFF_data,N2_N2_GAFF_data[2:3:end]);
plot!(dist_N2_N2_MMFF94S_data,N2_N2_MMFF94S_data[2:3:end]);
plot!(dist_N2_N2_CCSDT_data,N2_N2_CCSDT_data[2:3:end]);
plot!(dist_ecp_model,model_ecp_data);
plot!(ylims=(-4,2));
plot!(xlims=(2.5,5));
plot!(yticks=-4:2:2);
plot!(xticks=(2:5,[]));
plot!(ylabel=L"$\Delta E \ \left[ \mathrm{kJ/mol} \right]$");
plot!(legend=false);
annotate!(4.7,-3,text(L"\theta = \pi / 4",:center,:left,10));

fg2 = @pgf Axis(
    {
        legend_style =
        {
            at = Coordinate(0.5,1.5),
            anchor = "center",
            legend_columns = 2
        },
        xmajorgrids,
        ymajorgrids,
        width = 280, 
        height = 100,
        xmin = 2.75,
        xmax = 5,
        ymin = -4,
        ymax = 2,
        xtick = "3,3.5,4,4.5,5",
        xticklabels = "\\empty",
        ylabel = L"\Delta E \ \left[ \mathrm{kJ/mol} \right]",
        "grid style={line width=.1pt, draw=gray!10},major grid style={line width=.2pt,draw=gray!25}",
    },
    Plot(
        {
            "no marks",
            style = {"thick"},
            color = theme_palette(:auto)[1],
        },
        Coordinates(dist_ecp_model[model_ecp_data .< 10],model_ecp_data[model_ecp_data .< 10]),
        "node[anchor=center] at (4.75,-3) {\\fontsize{8pt}{8pt}\\selectfont \\color{black} \$ \\theta = \\pi / 4 \$}",
    ),
    Plot(
        {
            "no marks",
            style = {"thick"},
            color = theme_palette(:auto)[3],
        },
        Coordinates(dist_N2_N2_MMFF94S_data[N2_N2_MMFF94S_data[2:3:end] .< 5.0],N2_N2_MMFF94S_data[2:3:end][N2_N2_MMFF94S_data[2:3:end] .< 5.0]),
    ),
    Plot(
        {
            "no marks",
            style = {"thick"},
            color = theme_palette(:auto)[2],
        },
        Coordinates(dist_N2_N2_CCSDT_data,N2_N2_CCSDT_data[2:3:end]),
    ),
    Plot(
        {
            "no marks",
            style = {"thick"},
            color = theme_palette(:auto)[4],
        },
        Coordinates(dist_N2_N2_GAFF_data[N2_N2_GAFF_data[2:3:end] .< 5.0],N2_N2_GAFF_data[2:3:end][N2_N2_GAFF_data[2:3:end] .< 5.0]),
    )
);

# θ = π/2
model_ecp_data = zeros(Float64,0);
for dist in dist_ecp_model
    local mol_ecp_b = copy(mol_ecp_a);
    mol_len = abs(mol_ecp_b.atoms_data[1,1] - mol_ecp_b.atoms_data[2,1]);
    MoveAndRotateMolec!(mol_ecp_b,[0,0,0,dist/a0 + mol_len,0,0]);
    MoveAndRotateMolec!(mol_ecp_b,[0,0,π/2,0,0,0]);
    PolarizeMolecules!([mol_ecp_a,mol_ecp_b],xc_coeffs_ECP_Pol);

    local aux_energy = NaiveEnergyFromDensity([mol_ecp_a,mol_ecp_b]);
    local aux_energy += (XCEnergyFromDensity([mol_ecp_a,mol_ecp_b],order)*xc_coeffs_ECP)[1];
    push!(model_ecp_data,(aux_energy-model_ecp_e0)*kjmol*3.0);
end
N2_N2_ecp_model_data[3:3:end] = model_ecp_data;

p3 = plot(dist_N2_N2_GAFF_data,N2_N2_GAFF_data[3:3:end]);
plot!(dist_N2_N2_MMFF94S_data,N2_N2_MMFF94S_data[3:3:end]);
plot!(dist_N2_N2_CCSDT_data,N2_N2_CCSDT_data[3:3:end]);
plot!(dist_ecp_model,model_ecp_data);
plot!(ylims=(-4,2));
plot!(xlims=(2.5,5));
plot!(yticks=-4:2:2);
plot!(xlabel=L"$\Delta L \ \left[ \mathrm{\AA} \right]$");
plot!(ylabel=L"$\Delta E \ \left[ \mathrm{kJ/mol} \right]$");
plot!(legend=false);
annotate!(4.7,-3,text(L"\theta = \pi / 2",:center,:left,10));

s1 = 0.450;
s2 = 0.285;
s3 = 0.265;

N2_N2_P1 = plot(p1,p2,p3,layout=grid(3,1,heights=(s1,s2,s3)));
plot!(size=(550,410));
savefig("N2_N2_P1.pdf");

N2_N2_P2 = scatter(N2_N2_CCSDT_data,N2_N2_ecp_model_data,label="This Work (ECP Fit) (×3.0)")
scatter!(N2_N2_CCSDT_data,N2_N2_MMFF94S_data,label="MMFF94S (×3.0)");
scatter!(N2_N2_CCSDT_data,N2_N2_GAFF_data,label="GAFF");
plot!(xlims=(-50,1200),ylims=(-50,1200));
plot!(xticks=0:300:1200,yticks=0:300:1200);
plot!(legend=:bottomright);
plot!([-100,3000],[-100,3000],label=false);
annotate!(450,1100,text("N₂ + N₂",:center,:center,10));

F = @pgf Axis(
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
        height = 55,
    },
    Plot(
        {
            no_marks,
        },
        Coordinates([0],[0]),
        "node[] at (0.5,4.5) {\\includegraphics{/Users/joseantoniorodriguesromero/Documents/GitHub/PolarizedElectronDensityForceModel/N2+N2_diag.pdf}}"
    ),
);

FG = @pgf GroupPlot({group_style = { group_size = "1 by 4"}},F,fg1,fg2,fg3);

fg3 = @pgf Axis(
    {
        legend_style =
        {
            at = Coordinate(0.5,1.5),
            anchor = "center",
            legend_columns = 2
        },
        xmajorgrids,
        ymajorgrids,
        width = 280, 
        height = 100,
        xmin = 2.75,
        xmax = 5,
        ymin = -4,
        ymax = 2,
        xtick = "3,3.5,4,4.5,5",
        xticklabels = L"$3.0$,$3.5$,$4.0$,$4.5$,$5.0$",
        ylabel = L"\Delta E \ \left[ \mathrm{kJ/mol} \right]",
        xlabel = L"\Delta L \ \left[ \mathring{\mathrm{A}} \right]",
        "grid style={line width=.1pt, draw=gray!10},major grid style={line width=.2pt,draw=gray!25}",
    },
    Plot(
        {
            "no marks",
            style = {"thick"},
            color = theme_palette(:auto)[1],
        },
        Coordinates(dist_ecp_model[model_ecp_data .< 10],model_ecp_data[model_ecp_data .< 10]),
        "node[anchor=center] at (4.75,-3) {\\fontsize{8pt}{8pt}\\selectfont \\color{black} \$ \\theta = \\pi / 2 \$}",
    ),
    Plot(
        {
            "no marks",
            style = {"thick"},
            color = theme_palette(:auto)[3],
        },
        Coordinates(dist_N2_N2_MMFF94S_data[N2_N2_MMFF94S_data[3:3:end] .< 5.0],N2_N2_MMFF94S_data[3:3:end][N2_N2_MMFF94S_data[3:3:end] .< 5.0]),
    ),
    Plot(
        {
            "no marks",
            style = {"thick"},
            color = theme_palette(:auto)[2],
        },
        Coordinates(dist_N2_N2_CCSDT_data,N2_N2_CCSDT_data[3:3:end]),
    ),
    Plot(
        {
            "no marks",
            style = {"thick"},
            color = theme_palette(:auto)[4],
        },
        Coordinates(dist_N2_N2_GAFF_data[N2_N2_GAFF_data[3:3:end] .< 5.0],N2_N2_GAFF_data[3:3:end][N2_N2_GAFF_data[3:3:end] .< 5.0]),
    )
);

fp2 = @pgf Axis(
    {
        "legend pos = south east",
        xmajorgrids,
        ymajorgrids,
        width = 225, 
        height = 150,
        xmin = -60,
        xmax = 1200,
        ymin = -60,
        ymax = 1200,
        ytick = "0,300,600,900,1200",
        xtick = "0,300,600,900,1200",
        xticklabels = L"$0$,$300$,$600$,$900$,$1200$",
        yticklabels = L"$0$,$300$,$600$,$900$,$1200$",
        "grid style={line width=.1pt, draw=gray!10},major grid style={line width=.2pt,draw=gray!25}",
    },
    Plot(
        {
            "only marks",
            mark = "*",
            "mark options" = {"fill" = theme_palette(:auto)[1]},
        },
        Table([N2_N2_CCSDT_data[:],N2_N2_ecp_model_data[:]]),
        "node[] at (450,1050) {\$\\mathrm{N}_2 + \\mathrm{N}_2\$}",
    ),
    Plot(
        {
            "only marks",
            mark = "*",
            "mark options" = {"fill" = theme_palette(:auto)[2]},
        },
        Table([N2_N2_CCSDT_data[N2_N2_MMFF94S_data .< 1500],N2_N2_MMFF94S_data[N2_N2_MMFF94S_data .< 1500]]),
    ),
    Plot(
        {
            "only marks",
            mark = "*",
            "mark options" = {"fill" = theme_palette(:auto)[3]},
        },
        Table([N2_N2_CCSDT_data[N2_N2_GAFF_data .< 1500],N2_N2_GAFF_data[N2_N2_GAFF_data .< 1500]]),
    ),
    Plot(
        {
            no_markers,
            color = theme_palette(:auto)[4]
        },
        Coordinates([-200,4000],[-200,4000]),
    ),
    LegendEntry({anchor = "west"},L"\fontsize{8pt}{8pt}\selectfont $\mathrm{This \ Work} \ \left( \times 3.0 \right)$"),
    LegendEntry({anchor = "west"},L"\fontsize{8pt}{8pt}\selectfont $\mathrm{MMFF9S} \ \left( \times 3.0 \right)$"),
    LegendEntry({anchor = "west"},L"\fontsize{8pt}{8pt}\selectfont $\mathrm{GAFF}$"),
);