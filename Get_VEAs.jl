using Plots, SpecialFunctions, LaTeXStrings, Printf

include("GetPotentialFromAngles.jl")

function ReadXCCoeffs(file_name::String)
    # Reads the XC coefficients from the specified text file. The file is 
    # assumed to be in a single column format.
    ret = zeros(Float64,0);
    for line in readlines(file_name)
        push!(ret,parse(Float64,line));
    end

    return ret;
end

function ReadGaussianEnergy(file_name::String)
    for line in readlines(file_name)
        if contains(line,"CCSD(T)=")
            return parse(Float64,replace(split(line)[end],"D"=>"E"))
        end
    end

    return parse(Float64,"inf");
end

xc_coeffs_ECP = ReadXCCoeffs("ECP_XCCoeffs.txt");
xc_coeffs_FullE = ReadXCCoeffs("FullE_XCCoeffs.txt");
xc_coeffs_ECP_Pol = ReadXCCoeffs("ECP_XCCoeffs_Pol.txt");
xc_coeffs_FullE_Pol = ReadXCCoeffs("FullE_XCCoeffs_Pol.txt");

energies = zeros(Float64,4,3);

# H2 (Gaussian)
energies[1,1] += ReadGaussianEnergy("Test Data/Gaussian Data/Neutral/H2+H2/H2.log");
energies[1,1] -= ReadGaussianEnergy("Test Data/Gaussian Data/Anion/H2+H2/H2.log") / 2.0;

# H2 (ECP)
mol_a = ReadMolecule("H2_ecp_fitted_data.txt");
mol_b = ReadMolecule("H2_ecp_fitted_data.txt");

mol_a.charge = 0;
mol_b.charge = -1;

PolarizeMolecules!([mol_a],xc_coeffs_ECP_Pol);
PolarizeMolecules!([mol_b],xc_coeffs_ECP_Pol);

energies[1,2] -= NaiveEnergyFromDensity([mol_b]);
energies[1,2] += NaiveEnergyFromDensity([mol_a]);

energies[1,2] -= (XCEnergyFromDensity([mol_b],10)*xc_coeffs_ECP)[1];
energies[1,2] += (XCEnergyFromDensity([mol_a],10)*xc_coeffs_ECP)[1];

# H2 (Full E.)
mol_a = ReadMolecule("H2_fulle_fitted_data.txt");
mol_b = ReadMolecule("H2_fulle_fitted_data.txt");

mol_a.charge = 0;
mol_b.charge = -1;

PolarizeMolecules!([mol_a],xc_coeffs_FullE_Pol);
PolarizeMolecules!([mol_b],xc_coeffs_FullE_Pol);

energies[1,3] -= NaiveEnergyFromDensity([mol_b]);
energies[1,3] += NaiveEnergyFromDensity([mol_a]);

energies[1,3] -= (XCEnergyFromDensity([mol_b],10)*xc_coeffs_FullE)[1];
energies[1,3] += (XCEnergyFromDensity([mol_a],10)*xc_coeffs_FullE)[1];

# N2 (Gaussian)
energies[2,1] += ReadGaussianEnergy("Test Data/Gaussian Data/Neutral/N2+N2/N2.log");
energies[2,1] -= ReadGaussianEnergy("Test Data/Gaussian Data/Anion/N2+N2/N2.log");

# N2 (ECP)
mol_a = ReadMolecule("N2_ecp_fitted_data.txt");
mol_b = ReadMolecule("N2_ecp_fitted_data.txt");

mol_a.charge = 0;
mol_b.charge = -1;

PolarizeMolecules!([mol_a],xc_coeffs_ECP_Pol);
PolarizeMolecules!([mol_b],xc_coeffs_ECP_Pol);

energies[2,2] -= NaiveEnergyFromDensity([mol_b]);
energies[2,2] += NaiveEnergyFromDensity([mol_a]);

energies[2,2] -= (XCEnergyFromDensity([mol_b],10)*xc_coeffs_ECP)[1];
energies[2,2] += (XCEnergyFromDensity([mol_a],10)*xc_coeffs_ECP)[1];

# N2 (Full E.)
mol_a = ReadMolecule("N2_fulle_fitted_data.txt");
mol_b = ReadMolecule("N2_fulle_fitted_data.txt");

mol_a.charge = 0;
mol_b.charge = -1;

PolarizeMolecules!([mol_a],xc_coeffs_FullE_Pol);
PolarizeMolecules!([mol_b],xc_coeffs_FullE_Pol);

energies[2,3] -= NaiveEnergyFromDensity([mol_b]);
energies[2,3] += NaiveEnergyFromDensity([mol_a]);

energies[2,3] -= (XCEnergyFromDensity([mol_b],10)*xc_coeffs_FullE)[1];
energies[2,3] += (XCEnergyFromDensity([mol_a],10)*xc_coeffs_FullE)[1];

# O2 (Gaussian)
energies[3,1] += ReadGaussianEnergy("Test Data/Gaussian Data/Neutral/O2+O2 (2S+1=1)/O2.log");
energies[3,1] -= ReadGaussianEnergy("Test Data/Gaussian Data/Anion/O2+O2 (2S+1=2)/O2.log");

# O2 (ECP)
mol_a = ReadMolecule("O2_ecp_fitted_data.txt");
mol_b = ReadMolecule("O2_ecp_fitted_data.txt");

mol_a.charge = 0;
mol_b.charge = -1;

PolarizeMolecules!([mol_a],xc_coeffs_ECP_Pol);
PolarizeMolecules!([mol_b],xc_coeffs_ECP_Pol);

energies[3,2] -= NaiveEnergyFromDensity([mol_b]);
energies[3,2] += NaiveEnergyFromDensity([mol_a]);

energies[3,2] -= (XCEnergyFromDensity([mol_b],10)*xc_coeffs_ECP)[1];
energies[3,2] += (XCEnergyFromDensity([mol_a],10)*xc_coeffs_ECP)[1];

# O2 (Full E.)
mol_a = ReadMolecule("O2_fulle_fitted_data.txt");
mol_b = ReadMolecule("O2_fulle_fitted_data.txt");

mol_a.charge = 0;
mol_b.charge = -1;

PolarizeMolecules!([mol_a],xc_coeffs_FullE_Pol);
PolarizeMolecules!([mol_b],xc_coeffs_FullE_Pol);

energies[3,3] -= NaiveEnergyFromDensity([mol_b]);
energies[3,3] += NaiveEnergyFromDensity([mol_a]);

energies[3,3] -= (XCEnergyFromDensity([mol_b],10)*xc_coeffs_FullE)[1];
energies[3,3] += (XCEnergyFromDensity([mol_a],10)*xc_coeffs_FullE)[1];

# CO2 (Gaussian)
energies[4,1] += ReadGaussianEnergy("Test Data/Gaussian Data/Neutral/CO2+CO2/CO2.log");
energies[4,1] -= ReadGaussianEnergy("Test Data/Gaussian Data/Anion/CO2/CO2.log");

# CO2 (ECP)
mol_a = ReadMolecule("CO2_ecp_fitted_data.txt");
mol_b = ReadMolecule("CO2_ecp_fitted_data.txt");

mol_a.charge = 0;
mol_b.charge = -1;

PolarizeMolecules!([mol_a],xc_coeffs_ECP_Pol);
PolarizeMolecules!([mol_b],xc_coeffs_ECP_Pol);

energies[4,2] -= NaiveEnergyFromDensity([mol_b]);
energies[4,2] += NaiveEnergyFromDensity([mol_a]);

energies[4,2] -= (XCEnergyFromDensity([mol_b],10)*xc_coeffs_ECP)[1];
energies[4,2] += (XCEnergyFromDensity([mol_a],10)*xc_coeffs_ECP)[1];

# CO2 (FullE)
mol_a = ReadMolecule("CO2_fulle_fitted_data.txt");
mol_b = ReadMolecule("CO2_fulle_fitted_data.txt");

mol_a.charge = 0;
mol_b.charge = -1;

PolarizeMolecules!([mol_a],xc_coeffs_FullE_Pol);
PolarizeMolecules!([mol_b],xc_coeffs_FullE_Pol);

energies[4,3] -= NaiveEnergyFromDensity([mol_b]);
energies[4,3] += NaiveEnergyFromDensity([mol_a]);

energies[4,3] -= (XCEnergyFromDensity([mol_b],10)*xc_coeffs_FullE)[1];
energies[4,3] += (XCEnergyFromDensity([mol_a],10)*xc_coeffs_FullE)[1];

energies[:] .*= 27.211399;
energies = round.(energies,digits=3);