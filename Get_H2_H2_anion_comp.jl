using Plots, SpecialFunctions, LaTeXStrings
include("GetPotentialFromAngles.jl")

kjmol = 2625.5002;
a0 = 0.529177210903;
kCalMol_to_Hartree = 0.0015936011;

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

# Gaussian Data
function ReadGaussianEnergy(file_name::String)
    for line in readlines(file_name)
        if contains(line,"CCSD(T)=")
            return parse(Float64,replace(split(line)[end],"D"=>"E"))
        end
    end

    return parse(Float64,"inf");
end

dist_H2_H2_CCSDT_data = 1.0 .+ ((4.0/100.0).*collect(0:100));
H2_H2_CCSDT_data = zeros(Float64,303);
for i in 0:302
    local file_name = "Test Data/Gaussian Data/Anion/H2+H2/H2_H2_"*string(i)*".log"
    H2_H2_CCSDT_data[i+1] = ReadGaussianEnergy(file_name)
end

# θ = 0
H2_H2_ecp_model_data = zeros(Float64,303);
H2_H2_fulle_model_data = zeros(Float64,303);
dist_ecp_model = (1.0 .+ ((4.0/100.0).*collect(0:100)));
model_ecp_data = zeros(Float64,0);
mol_ecp_a = ReadMolecule("H2_ecp_fitted_data.txt");

PolarizeMolecules!([mol_ecp_a],xc_coeffs_ECP_Pol);
order = ceil(Int,length(xc_coeffs_ECP) / 4) - 1;

for dist in dist_ecp_model
    local mol_ecp_b = copy(mol_ecp_a);
    mol_ecp_b.charge = -1;
    mol_len = abs(mol_ecp_b.atoms_data[1,1] - mol_ecp_b.atoms_data[2,1]);
    MoveAndRotateMolec!(mol_ecp_b,[0,0,0,dist/a0 + mol_len,0,0]);
    PolarizeMolecules!([mol_ecp_a,mol_ecp_b],xc_coeffs_ECP_Pol);

    local aux_energy = NaiveEnergyFromDensity([mol_ecp_a,mol_ecp_b]);
    local aux_energy += (XCEnergyFromDensity([mol_ecp_a,mol_ecp_b],order)*xc_coeffs_ECP)[1];
    push!(model_ecp_data,(aux_energy));
end
H2_H2_ecp_model_data[1:3:end] = model_ecp_data;

# θ = π/4
model_ecp_data = zeros(Float64,0);
for dist in dist_ecp_model
    local mol_ecp_b = copy(mol_ecp_a);
    mol_ecp_b.charge = -1;
    mol_len = abs(mol_ecp_b.atoms_data[1,1] - mol_ecp_b.atoms_data[2,1]);
    MoveAndRotateMolec!(mol_ecp_b,[0,0,0,dist/a0 + mol_len,0,0]);
    MoveAndRotateMolec!(mol_ecp_b,[0,0,π/4,0,0,0]);
    PolarizeMolecules!([mol_ecp_a,mol_ecp_b],xc_coeffs_ECP_Pol);

    local aux_energy = NaiveEnergyFromDensity([mol_ecp_a,mol_ecp_b]);
    local aux_energy += (XCEnergyFromDensity([mol_ecp_a,mol_ecp_b],order)*xc_coeffs_ECP)[1];
    push!(model_ecp_data,(aux_energy));
end
H2_H2_ecp_model_data[2:3:end] = model_ecp_data;

# θ = π/2
model_ecp_data = zeros(Float64,0);
for dist in dist_ecp_model
    local mol_ecp_b = copy(mol_ecp_a);
    mol_ecp_b.charge = -1;
    mol_len = abs(mol_ecp_b.atoms_data[1,1] - mol_ecp_b.atoms_data[2,1]);
    MoveAndRotateMolec!(mol_ecp_b,[0,0,0,dist/a0 + mol_len,0,0]);
    MoveAndRotateMolec!(mol_ecp_b,[0,0,π/2,0,0,0]);
    PolarizeMolecules!([mol_ecp_a,mol_ecp_b],xc_coeffs_ECP_Pol);

    local aux_energy = NaiveEnergyFromDensity([mol_ecp_a,mol_ecp_b]);
    local aux_energy += (XCEnergyFromDensity([mol_ecp_a,mol_ecp_b],order)*xc_coeffs_ECP)[1];
    push!(model_ecp_data,(aux_energy));
end
H2_H2_ecp_model_data[3:3:end] = model_ecp_data;

H2_H2_ecp_model_data .-= H2_H2_ecp_model_data[1:3:end][end];
H2_H2_CCSDT_data .-= H2_H2_CCSDT_data[1:3:end][end];

H2_H2_ecp_model_data .*= kjmol;
H2_H2_CCSDT_data .*= kjmol;

H2_H2_ecp_model_data = H2_H2_ecp_model_data[isinf.(H2_H2_CCSDT_data) .== false];
H2_H2_CCSDT_data = H2_H2_CCSDT_data[isinf.(H2_H2_CCSDT_data) .== false];

M = hcat([minimum(H2_H2_CCSDT_data),1.0],[maximum(H2_H2_CCSDT_data),1.0])';
Y = [minimum(H2_H2_ecp_model_data),maximum(H2_H2_ecp_model_data)];
X = M \ Y;

if X[2] >= 0
    aux_label = "f(x) = "*string(round(X[1],digits=4))*"*x + "*string(round(X[2],digits=4));
else
    aux_label = "f(x) = "*string(round(X[1],digits=4))*"*x − "*string(round(-X[2],digits=4));
end

X0 = -55;
X1 = 150;

Y0 = -10;
Y1 = 350;

H2_H2_P2 = scatter(H2_H2_CCSDT_data,H2_H2_ecp_model_data,label=false);
plot!(xlims=(X0,X1),ylims=(Y0,Y1));
plot!(xticks=-50:50:150,yticks=0:70:350);
plot!([X0,X1],[X0*X[1]+X[2],X1*X[1]+X[2]],label=aux_label);
plot!(ylabel="FF Model\n[kJ/mol]");
plot!(legend=:topleft);
plot!(yguidefontsize=8);

px = X0 + 0.875*(X1-X0);
py = Y0 + 0.1*(Y1-Y0);
annotate!(px,py,text("(H₂ + H₂)⁻",:center,:center,10));

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

aux_legend = "\\fontsize{8pt}{8pt}\\selectfont \$f \\left( x \\right) = ";
aux_legend = aux_legend * string(round(X[1],digits=3)) * "x + ";
aux_legend = aux_legend * string(round(X[2],digits=3)) * "\$"

fp1 = @pgf Axis(
    {
        "legend pos = north west",
        xmajorgrids,
        ymajorgrids,
        width = 225, 
        height = 150,
        xmin = -70,
        xmax = 150,
        ymin = -30,
        ymax = 350,
        xtick = "-50,0,50,100,150",
        ytick = "0,70,140,210,280,350",
        ylabel = L"$\mathrm{FF \ Model} \ \left[ \mathrm{kJ/mol} \right]$",
        "grid style={line width=.1pt, draw=gray!10},major grid style={line width=.2pt,draw=gray!25}",
    },
    Plot(
        {
            "only marks",
            mark = "*",
            "mark options" = {"fill" = theme_palette(:auto)[1]},
            "forget plot",
        },
        Table([H2_H2_CCSDT_data[:],H2_H2_ecp_model_data[:]]),
        "node[] at (110,30) {\$ \\left( \\mathrm{H}_2 + \\mathrm{H}_2 \\right)^- \$}",
    ),
    Plot(
        {
            no_markers,
            color = theme_palette(:auto)[2]
        },
        Coordinates([-200,4000],[X[1]*(-200)+X[2],X[1]*(4000)+X[2]]),
    ),
    LegendEntry({anchor = "west"},aux_legend),
)