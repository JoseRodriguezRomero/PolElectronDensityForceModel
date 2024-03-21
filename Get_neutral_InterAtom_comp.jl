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

# Auxiliary Functions
function ReadGaussianEnergy(file_name::String)
    for line in readlines(file_name)
        if contains(line,"CCSD(T)=")
            return parse(Float64,replace(split(line)[end],"D"=>"E"))
        end
    end

    return parse(Float64,"inf");
end

function ECP_AtomFromAtomicNum(atomic_num::Int)
    all_coeffs = readlines("AtomsCoeffs_ECP.txt");
    coeff_inf = split(all_coeffs[atomic_num+1]);

    atom = Molecule(Float64);
    atoms_data = zeros(Float64,1,4);
    atoms_data[end] = parse(Float64,coeff_inf[3]);
    atom.atoms_data = atoms_data;

    cloud_data = zeros(Float64,3,6);
    cloud_data[:,4] = parse.(Float64,coeff_inf[4:2:end]);
    cloud_data[:,5] = parse.(Float64,coeff_inf[5:2:end]);
    cloud_data[:,6] .= 1.0;
    atom.cloud_data = cloud_data;

    return atom;
end

function FullE_AtomFromAtomicNum(atomic_num::Int)
    all_coeffs = readlines("AtomsCoeffs_FullE.txt");
    coeff_inf_1 = split(all_coeffs[3*(atomic_num-1)+1+1]);
    coeff_inf_2 = split(all_coeffs[3*(atomic_num-1)+2+1]);
    coeff_inf_3 = split(all_coeffs[3*(atomic_num-1)+3+1]);

    atom = Molecule(Float64);
    atoms_data = zeros(Float64,1,4);
    atoms_data[end] = parse(Float64,coeff_inf_1[3]);
    atom.atoms_data = atoms_data;

    cloud_data = zeros(Float64,9,6);
    c_coeffs = parse.(Float64,coeff_inf_1[4:2:end]);
    c_coeffs = vcat(c_coeffs,parse.(Float64,coeff_inf_2[1:2:end]));
    c_coeffs = vcat(c_coeffs,parse.(Float64,coeff_inf_3[1:2:end]));

    λ_coeffs = parse.(Float64,coeff_inf_1[5:2:end]);
    λ_coeffs = vcat(λ_coeffs,parse.(Float64,coeff_inf_2[2:2:end]));
    λ_coeffs = vcat(λ_coeffs,parse.(Float64,coeff_inf_3[2:2:end]));

    cloud_data[:,4] = c_coeffs;
    cloud_data[:,5] = λ_coeffs;
    cloud_data[:,6] .= 1.0;
    atom.cloud_data = cloud_data;

    return atom;
end

# H2 Gaussian Data
H2_CCSDT_dist = 0.5 .+ (3.0/100.0).*(0:100);
H2_CCSDT_energy = zeros(Float64,101);

H_CCSDT_energy = ReadGaussianEnergy("Test Data/Gaussian Data/Neutral/H2+H2/H2.log");

for i in 0:100
    local file_name = "Test Data/Gaussian Data/Neutral/H2/H2_"*string(i)*".log";
    H2_CCSDT_energy[i+1] = ReadGaussianEnergy(file_name) - H_CCSDT_energy;
end

# H2 Model Data
H2_model_dist = 0.5 .+ (3.0/100.0).*(0:100);
H2_model_energy = zeros(Float64,101);

H1 = ECP_AtomFromAtomicNum(1);
H2 = ECP_AtomFromAtomicNum(1);

MoveAndRotateMolec!(H2,[0,0,0,0.74/a0,0,0]);
PolarizeMolecules!([H1,H2],xc_coeffs_ECP_Pol);

H2_0_model_energy = NaiveEnergyFromDensity([H1,H2]);
H2_0_model_energy += (XCEnergyFromDensity([H1,H2],10)*xc_coeffs_ECP)[1];

for i in 0:100
    local H1 = ECP_AtomFromAtomicNum(1);
    local H2 = ECP_AtomFromAtomicNum(1);

    dx = (0.5 .+ (3.0/100.0)*i)/a0;

    MoveAndRotateMolec!(H2,[0,0,0,dx,0,0]);
    PolarizeMolecules!([H1,H2],xc_coeffs_ECP_Pol);
    
    local aux_energy = NaiveEnergyFromDensity([H1,H2]);
    local aux_energy += (XCEnergyFromDensity([H1,H2],10)*xc_coeffs_ECP)[1];
    H2_model_energy[i+1] = aux_energy - H2_0_model_energy;
end

# N2 Gaussian Data
N2_CCSDT_dist = 0.8 .+ (3.0/100.0).*(0:100);
N2_CCSDT_energy = zeros(Float64,101);

N2_0_CCSDT_energy = ReadGaussianEnergy("Test Data/Gaussian Data/Neutral/N2+N2/N2.log");

for i in 0:100
    local file_name = "Test Data/Gaussian Data/Neutral/N2/N2_"*string(i)*".log";
    N2_CCSDT_energy[i+1] = ReadGaussianEnergy(file_name) - N2_0_CCSDT_energy;
end

N2_CCSDT_dist = N2_CCSDT_dist[isinf.(N2_CCSDT_energy) .== false];
N2_CCSDT_energy = N2_CCSDT_energy[isinf.(N2_CCSDT_energy) .== false];

# N2 Model Data
N2_ecp_model_dist = 0.8 .+ (3.0/100.0).*(0:100);
N2_ecp_model_energy = zeros(Float64,101);

N1 = ECP_AtomFromAtomicNum(7);
N2 = ECP_AtomFromAtomicNum(7);

MoveAndRotateMolec!(N2,[0,0,0,1.09/a0,0,0]);
PolarizeMolecules!([N1,N2],xc_coeffs_ECP_Pol);

N2_0_model_energy = NaiveEnergyFromDensity([N1,N2]);
N2_0_model_energy += (XCEnergyFromDensity([N1,N2],10)*xc_coeffs_ECP)[1];

for i in 0:100
    local N1 = ECP_AtomFromAtomicNum(7);
    local N2 = ECP_AtomFromAtomicNum(7);

    dx = (0.8 .+ (3.0/100.0)*i)/a0;

    MoveAndRotateMolec!(N2,[0,0,0,dx,0,0]);
    PolarizeMolecules!([N1,N2],xc_coeffs_ECP_Pol);
    
    local aux_energy = NaiveEnergyFromDensity([N1,N2]);
    local aux_energy += (XCEnergyFromDensity([N1,N2],10)*xc_coeffs_ECP)[1];
    N2_ecp_model_energy[i+1] = aux_energy - N2_0_model_energy;
end

N2_fulle_model_dist = 0.8 .+ (3.0/100.0).*(0:100);
N2_fulle_model_energy = zeros(Float64,101);

N1 = FullE_AtomFromAtomicNum(7);
N2 = FullE_AtomFromAtomicNum(7);

MoveAndRotateMolec!(N2,[0,0,0,1.09/a0,0,0]);
PolarizeMolecules!([N1,N2],xc_coeffs_FullE_Pol);

N2_0_model_energy = NaiveEnergyFromDensity([N1,N2]);
N2_0_model_energy += (XCEnergyFromDensity([N1,N2],10)*xc_coeffs_FullE)[1];

for i in 0:100
    local N1 = FullE_AtomFromAtomicNum(7);
    local N2 = FullE_AtomFromAtomicNum(7);

    dx = (0.8 .+ (3.0/100.0)*i)/a0;

    MoveAndRotateMolec!(N2,[0,0,0,dx,0,0]);
    PolarizeMolecules!([N1,N2],xc_coeffs_FullE_Pol);
    
    local aux_energy = NaiveEnergyFromDensity([N1,N2]);
    local aux_energy += (XCEnergyFromDensity([N1,N2],10)*xc_coeffs_FullE)[1];
    N2_fulle_model_energy[i+1] = aux_energy - N2_0_model_energy;
end

# O2 Gaussian Data
O2_CCSDT_dist = 0.8 .+ (3.0/100.0).*(0:100);
O2_CCSDT_energy = zeros(Float64,101);

# O_CCSDT_energy = ReadGaussianEnergy("Test Data/Gaussian Data/Neutral/O2+O2 (2S+1=3)/O2.log");

for i in 0:100
    local file_name = "Test Data/Gaussian Data/Neutral/O2/O2_"*string(i)*".log";
    # O2_CCSDT_energy[i+1] = ReadGaussianEnergy(file_name) - O_CCSDT_energy;
    O2_CCSDT_energy[i+1] = ReadGaussianEnergy(file_name);
end

O2_CCSDT_energy = O2_CCSDT_energy .- minimum(O2_CCSDT_energy);

O2_CCSDT_dist = O2_CCSDT_dist[isinf.(O2_CCSDT_energy) .== false];
O2_CCSDT_energy = O2_CCSDT_energy[isinf.(O2_CCSDT_energy) .== false];

# O2 Model Data
O2_ecp_model_dist = 0.8 .+ (3.0/100.0).*(0:100);
O2_ecp_model_energy = zeros(Float64,101);

O1 = ECP_AtomFromAtomicNum(8);
O2 = ECP_AtomFromAtomicNum(8);

MoveAndRotateMolec!(O2,[0,0,0,1.09/a0,0,0]);
PolarizeMolecules!([O1,O2],xc_coeffs_ECP_Pol);

O2_0_model_energy = NaiveEnergyFromDensity([O1,O2]);
O2_0_model_energy += (XCEnergyFromDensity([O1,O2],10)*xc_coeffs_ECP)[1];

for i in 0:100
    local O1 = ECP_AtomFromAtomicNum(8);
    local O2 = ECP_AtomFromAtomicNum(8);

    dx = (0.8 .+ (3.0/100.0)*i)/a0;

    MoveAndRotateMolec!(O2,[0,0,0,dx,0,0]);
    PolarizeMolecules!([O1,O2],xc_coeffs_ECP_Pol);
    
    local aux_energy = NaiveEnergyFromDensity([O1,O2]);
    local aux_energy += (XCEnergyFromDensity([O1,O2],10)*xc_coeffs_ECP)[1];
    O2_ecp_model_energy[i+1] = aux_energy - O2_0_model_energy;
end

O2_fulle_model_dist = 0.8 .+ (3.0/100.0).*(0:100);
O2_fulle_model_energy = zeros(Float64,101);

O1 = FullE_AtomFromAtomicNum(8);
O2 = FullE_AtomFromAtomicNum(8);

MoveAndRotateMolec!(O2,[0,0,0,1.09/a0,0,0]);
PolarizeMolecules!([O1,O2],xc_coeffs_FullE_Pol);

O2_0_model_energy = NaiveEnergyFromDensity([O1,O2]);
O2_0_model_energy += (XCEnergyFromDensity([O1,O2],10)*xc_coeffs_FullE)[1];

for i in 0:100
    local O1 = FullE_AtomFromAtomicNum(8);
    local O2 = FullE_AtomFromAtomicNum(8);

    dx = (0.8 .+ (3.0/100.0)*i)/a0;

    MoveAndRotateMolec!(O2,[0,0,0,dx,0,0]);
    PolarizeMolecules!([O1,O2],xc_coeffs_FullE_Pol);
    
    local aux_energy = NaiveEnergyFromDensity([O1,O2]);
    local aux_energy += (XCEnergyFromDensity([O1,O2],10)*xc_coeffs_FullE)[1];
    O2_fulle_model_energy[i+1] = aux_energy - O2_0_model_energy;
end

# B + Ar Gaussian Data
B_Ar_CCSDT_dist = 0.8 .+ (3.0/100.0).*(0:100);
B_Ar_CCSDT_energy = zeros(Float64,101);

B_CCSDT_energy = ReadGaussianEnergy("Test Data/Gaussian Data/Neutral/B+Ar/B.log");
Ar_CCSDT_energy = ReadGaussianEnergy("Test Data/Gaussian Data/Neutral/B+Ar/Ar.log");

for i in 0:100
    local file_name = "Test Data/Gaussian Data/Neutral/B+Ar/B_Ar_"*string(i)*".log";
    B_Ar_CCSDT_energy[i+1] = ReadGaussianEnergy(file_name) - B_CCSDT_energy - Ar_CCSDT_energy;
end

B_Ar_CCSDT_dist = B_Ar_CCSDT_dist[isinf.(B_Ar_CCSDT_energy) .== false];
B_Ar_CCSDT_energy = B_Ar_CCSDT_energy[isinf.(B_Ar_CCSDT_energy) .== false];

# B + Ar Model Data
B_Ar_ecp_model_dist = 0.8 .+ (3.0/100.0).*(0:100);
B_Ar_ecp_model_energy = zeros(Float64,101);

B1 = ECP_AtomFromAtomicNum(5);
B_ecp_model_energy = NaiveEnergyFromDensity([B1]);
B_ecp_model_energy += (XCEnergyFromDensity([B1],10)*xc_coeffs_ECP)[1];

Ar2 = ECP_AtomFromAtomicNum(18);
Ar_ecp_model_energy = NaiveEnergyFromDensity([Ar2]);
Ar_ecp_model_energy += (XCEnergyFromDensity([Ar2],10)*xc_coeffs_ECP)[1];

for i in 0:100
    local B1 = ECP_AtomFromAtomicNum(5);
    local Ar2 = ECP_AtomFromAtomicNum(18);

    dx = (0.8 .+ (3.0/100.0)*i)/a0;

    MoveAndRotateMolec!(Ar2,[0,0,0,dx,0,0]);
    PolarizeMolecules!([B1,Ar2],xc_coeffs_ECP_Pol);
    
    local aux_energy = NaiveEnergyFromDensity([B1,Ar2]);
    local aux_energy += (XCEnergyFromDensity([B1,Ar2],10)*xc_coeffs_ECP)[1];
    B_Ar_ecp_model_energy[i+1] = aux_energy - B_ecp_model_energy - Ar_ecp_model_energy;
end

B_Ar_fulle_model_dist = 0.8 .+ (3.0/100.0).*(0:100);
B_Ar_fulle_model_energy = zeros(Float64,101);

B1 = FullE_AtomFromAtomicNum(5);
B_fulle_model_energy = NaiveEnergyFromDensity([B1]);
B_fulle_model_energy += (XCEnergyFromDensity([B1],10)*xc_coeffs_FullE)[1];

Ar2 = FullE_AtomFromAtomicNum(18);
Ar_fulle_model_energy = NaiveEnergyFromDensity([Ar2]);
Ar_fulle_model_energy += (XCEnergyFromDensity([Ar2],10)*xc_coeffs_FullE)[1];

for i in 0:100
    global B1 = FullE_AtomFromAtomicNum(5);
    global Ar2 = FullE_AtomFromAtomicNum(18);

    dx = (0.8 .+ (3.0/100.0)*i)/a0;

    MoveAndRotateMolec!(Ar2,[0,0,0,dx,0,0]);
    PolarizeMolecules!([B1,Ar2],xc_coeffs_FullE_Pol);
    
    local aux_energy = NaiveEnergyFromDensity([B1,Ar2]);
    local aux_energy += (XCEnergyFromDensity([B1,Ar2],10)*xc_coeffs_FullE)[1];
    B_Ar_fulle_model_energy[i+1] = aux_energy - B_fulle_model_energy - Ar_fulle_model_energy;
end

H2_CCSDT_energy .*= kjmol;
H2_model_energy .*= kjmol;

N2_CCSDT_energy .*= kjmol;
N2_ecp_model_energy .*= kjmol;
N2_fulle_model_energy .*= kjmol;

O2_CCSDT_energy .*= kjmol;
O2_ecp_model_energy .*= kjmol;
O2_fulle_model_energy .*= kjmol;

B_Ar_CCSDT_energy .*= kjmol;
B_Ar_ecp_model_energy .*= kjmol;
B_Ar_fulle_model_energy .*= kjmol;

# Plot all the data
p1 = plot(H2_CCSDT_dist,H2_CCSDT_energy,label="CCSD(T)/cc-pVTZ");
plot!(H2_model_dist,H2_model_energy,label="This Work (ECP Fit)");
plot!([1,2],[10000,10000],label="This Work (Full E. Fit)");
plot!(xlims=(0.5,3.5),ylims=(-600,600));
plot!(xticks=(0.5:0.5:3.5,[]),yticks=(-600:400:600));
plot!(legend=:outertop,legendcolumns=3);
plot!(ylabel=L"$\Delta E \quad \left[ \mathrm{kJ/mol} \right]$");

fg1 = @pgf Axis(
    {
        legend_style =
        {
            at = Coordinate(0.5,1.4),
            anchor = "center",
            legend_columns = 3
        },
        xmajorgrids,
        ymajorgrids,
        width = 280, 
        height = 100,
        xmin = 0.5,
        xmax = 3.5,
        ymin = -600,
        ymax = 600,
        xtick = "0.5,1,1.5,2,2.5,3",
        ytick = "-600,-200,200,600",
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
        Coordinates(H2_CCSDT_dist,H2_CCSDT_energy),
        "node[anchor=center] at (3.25,-400) {\\fontsize{8pt}{8pt}\\selectfont \\color{black} \$ \\mathrm{H} + \\mathrm{H} \$}",
    ),
    Plot(
        {
            "no marks",
            style = {"thick"},
            color = theme_palette(:auto)[2],
        },
        Coordinates(H2_model_dist,H2_model_energy),
    ),
    Plot(
        {
            "no marks",
            style = {"thick"},
            color = theme_palette(:auto)[3],
        },
        Coordinates([0],[-700]),
    ),
    LegendEntry({anchor = "west"},L"\fontsize{8pt}{8pt}\selectfont $\mathrm{CCSD(T)/cc\text{-}pVTZ}$"),
    LegendEntry({anchor = "west"},L"\fontsize{8pt}{8pt}\selectfont $\mathrm{This \ Work \ (ECP \ fit)}$"),
    LegendEntry({anchor = "west"},L"\fontsize{8pt}{8pt}\selectfont $\mathrm{This \ Work \ (Full E. \ fit)}$"),
);

px = 3.25;
py = -600 + 0.15*(600-(-600));
annotate!(px,py,text("H + H",:center,:center,8));

p2 = plot(N2_CCSDT_dist,N2_CCSDT_energy);
plot!(N2_ecp_model_dist,N2_ecp_model_energy);
plot!(N2_fulle_model_dist,N2_fulle_model_energy);
plot!(xlims=(0.5,3.5),ylims=(-2000,2000));
plot!(xticks=(0.5:0.5:3.5,[]),yticks=-2000:1000:2000);
plot!(legend=false);
plot!(ylabel=L"$\Delta E \quad \left[ \mathrm{kJ/mol} \right]$");

fg2 = @pgf Axis(
    {
        legend_style =
        {
            at = Coordinate(0.5,1.4),
            anchor = "center",
            legend_columns = 3
        },
        xmajorgrids,
        ymajorgrids,
        width = 280, 
        height = 100,
        xmin = 0.5,
        xmax = 3.5,
        ymin = -2000,
        ymax = 2000,
        xtick = "0.5,1,1.5,2,2.5,3",
        xticklabels = "\\empty",
        ytick = "-2000,-1000,0,1000,2000",
        yticklabels = L"$-2000$,$-1000$,$0$,$1000$,$2000$",
        ylabel = L"\Delta E \ \left[ \mathrm{kJ/mol} \right]",
        "grid style={line width=.1pt, draw=gray!10},major grid style={line width=.2pt,draw=gray!25}",
    },
    Plot(
        {
            "no marks",
            style = {"thick"},
            color = theme_palette(:auto)[1],
        },
        Coordinates(N2_CCSDT_dist,N2_CCSDT_energy),
        "node[anchor=center] at (3.25,1500) {\\fontsize{8pt}{8pt}\\selectfont \\color{black} \$ \\mathrm{N} + \\mathrm{N} \$}",
    ),
    Plot(
        {
            "no marks",
            style = {"thick"},
            color = theme_palette(:auto)[2],
        },
        Coordinates(N2_ecp_model_dist,N2_ecp_model_energy),
    ),
    Plot(
        {
            "no marks",
            style = {"thick"},
            color = theme_palette(:auto)[3],
        },
        Coordinates(N2_fulle_model_dist,N2_fulle_model_energy),
    ),
);

px = 3.25;
py = -2000 + 0.875*(2000-(-2000));
annotate!(px,py,text("N + N",:center,:center,8));

p3 = plot(O2_CCSDT_dist,O2_CCSDT_energy);
plot!(O2_ecp_model_dist,O2_ecp_model_energy);
plot!(O2_fulle_model_dist,O2_fulle_model_energy);
plot!(xlims=(0.5,3.5),ylims=(-4000,4000));
plot!(xticks=(0.5:0.5:3.5,[]),yticks=-4000:2000:4000);
plot!(legend=false);
plot!(ylabel=L"$\Delta E \quad \left[ \mathrm{kJ/mol} \right]$");

fg3 = @pgf Axis(
    {
        legend_style =
        {
            at = Coordinate(0.5,1.4),
            anchor = "center",
            legend_columns = 3
        },
        xmajorgrids,
        ymajorgrids,
        width = 280, 
        height = 100,
        xmin = 0.5,
        xmax = 3.5,
        ymin = -4000,
        ymax = 4000,
        xtick = "0.5,1,1.5,2,2.5,3",
        xticklabels = "\\empty",
        ytick = "-4000,-2000,0,2000,4000",
        yticklabels = L"$-4000$,$-2000$,$0$,$2000$,$4000$",
        ylabel = L"\Delta E \ \left[ \mathrm{kJ/mol} \right]",
        "grid style={line width=.1pt, draw=gray!10},major grid style={line width=.2pt,draw=gray!25}",
    },
    Plot(
        {
            "no marks",
            style = {"thick"},
            color = theme_palette(:auto)[1],
        },
        Coordinates(O2_CCSDT_dist,O2_CCSDT_energy),
        "node[anchor=center] at (3.25,3000) {\\fontsize{8pt}{8pt}\\selectfont \\color{black} \$ \\mathrm{O} + \\mathrm{O} \$}",
    ),
    Plot(
        {
            "no marks",
            style = {"thick"},
            color = theme_palette(:auto)[2],
        },
        Coordinates(O2_ecp_model_dist,O2_ecp_model_energy),
    ),
    Plot(
        {
            "no marks",
            style = {"thick"},
            color = theme_palette(:auto)[3],
        },
        Coordinates(O2_fulle_model_dist,O2_fulle_model_energy),
    ),
);

px = 3.25;
py = -4000 + 0.875*(4000-(-4000));
annotate!(px,py,text("O + O",:center,:center,8));

p4 = plot(B_Ar_CCSDT_dist,B_Ar_CCSDT_energy);
plot!(B_Ar_ecp_model_dist,B_Ar_ecp_model_energy);
plot!(B_Ar_fulle_model_dist,B_Ar_fulle_model_energy);
plot!(xlims=(0.5,3.5),ylims=(1.0E-2,1.0E10),yaxis=:log);
plot!(xticks=(0.5:0.5:3.5),yticks=(10.0.^collect(-2:4:10)));
plot!(legend=false,xlabel=L"$\Delta L \quad [\AA]$");
plot!(ylabel=L"$\Delta E \quad \left[ \mathrm{kJ/mol} \right]$");

fg4 = @pgf Axis(
    {
        legend_style =
        {
            at = Coordinate(0.5,1.4),
            anchor = "center",
            legend_columns = 3
        },
        xmajorgrids,
        ymajorgrids,
        "ymode = log",
        width = 280, 
        height = 100,
        xmin = 0.5,
        xmax = 3.5,
        ymin = 1.0E-2,
        ymax = 1.0E10,
        xtick = "0.5,1,1.5,2,2.5,3,3.5",
        xticklabels = L"$0.5$,$1.0$,$1.5$,$2.0$,$2.5$,$3.0$,$3.5$",
        ytick = "1.09E-2,1.0E2,1.0E6,1.0E10",
        yticklabels = L"$10^{-2}$,$10^{2}$,$10^{6}$,$10^{10}$",
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
        Coordinates(B_Ar_CCSDT_dist,B_Ar_CCSDT_energy),
        "node[anchor=center] at (3.25,1.0E4) {\\fontsize{8pt}{8pt}\\selectfont \\color{black} \$ \\mathrm{B} + \\mathrm{Ar} \$}",
    ),
    Plot(
        {
            "no marks",
            style = {"thick"},
            color = theme_palette(:auto)[2],
        },
        Coordinates(B_Ar_ecp_model_dist,B_Ar_ecp_model_energy),
    ),
    Plot(
        {
            "no marks",
            style = {"thick"},
            color = theme_palette(:auto)[3],
        },
        Coordinates(B_Ar_fulle_model_dist,B_Ar_fulle_model_energy),
    ),
);

px = 3.25;
py = -750 + 0.85*(6000-(-750));
annotate!(px,py,text("B + Ar",:center,:center,8));

pAll = plot(p1,p2,p3,p4,layout=grid(4,1,heights=(0.31,0.24,0.24,0.23)));
plot!(size=(600,550))

FG = @pgf GroupPlot({group_style = { group_size = "1 by 4"}},fg1,fg2,fg3,fg4);