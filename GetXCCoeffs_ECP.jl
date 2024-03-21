using Printf, Distributed, BenchmarkTools

if nworkers() != 11
    addprocs(12-nprocs());
end

@everywhere include("GetPotentialFromAngles.jl")
@everywhere include("ProcessFittingData.jl")

@everywhere H2O_gauss_e0 = ReadGaussianFile("Training Data/Gaussian Data/ECP/H2O.log","CCSD(T)");

@everywhere all_H2O_H2O_gauss_e = zeros(Float64,10000,1);
@everywhere all_H2O_H2O_angs_and_disps = zeros(Float64,10000,6);

# OH Coord
@everywhere all_angs_and_disps_data = readlines("Training Data/Gaussian Data/ECP/OH Coord/RotationsAndDisplacements.txt");
@everywhere for i in 1:5000
    gauss_file = "Training Data/Gaussian Data/ECP/OH Coord/H2O_H2O_"
    gauss_file *= string(i-1)*".log";

    all_H2O_H2O_gauss_e[i] = ReadGaussianFile(gauss_file,"CCSD(T)");

    angs_and_disps_data = split(all_angs_and_disps_data[i+1]);
    for j in 1:6
        all_H2O_H2O_angs_and_disps[i,j] = parse(Float64,angs_and_disps_data[j]);
    end
end

# OO Coord
@everywhere all_angs_and_disps_data = readlines("Training Data/Gaussian Data/ECP/OO Coord/RotationsAndDisplacements.txt");
@everywhere for i in 1:5000
    gauss_file = "Training Data/Gaussian Data/ECP/OO Coord/H2O_H2O_"
    gauss_file *= string(i-1)*".log";

    all_H2O_H2O_gauss_e[i+5000] = ReadGaussianFile(gauss_file,"CCSD(T)");

    angs_and_disps_data = split(all_angs_and_disps_data[i+1]);
    for j in 1:6
        all_H2O_H2O_angs_and_disps[i+5000,j] = parse(Float64,angs_and_disps_data[j]);
    end
end

@everywhere ecp_xc_pol_coeffs = zeros(Float64,0);
@everywhere for coeff in readlines("ECP_XCCoeffs_Pol.txt");
    push!(ecp_xc_pol_coeffs,parse(Float64,coeff));
end

@everywhere order = 10;
@everywhere mol_a = ReadMolecule("H2O_ECP_fitted_data.txt");
@everywhere PolarizeMolecules!([mol_a],ecp_xc_pol_coeffs);
@everywhere H2O_e0_naive = NaiveEnergyFromDensity([mol_a]);
@everywhere H2O_e0_XC = XCEnergyFromDensity([mol_a],order);

@everywhere function Foo(thread_id::Int)
    auxM = zeros(Float64,4,4);
    auxY = zeros(Float64,4,1);

    for ii in thread_id:nworkers():10000
        mol_b = copy(mol_a);
        angs_and_disps = all_H2O_H2O_angs_and_disps[ii,:];
        MoveAndRotateMolec!(mol_b,angs_and_disps);
    
        molecs = [mol_a,mol_b];
        PolarizeMolecules!(molecs,ecp_xc_pol_coeffs);
    
        gauss_e = all_H2O_H2O_gauss_e[ii];
        naive_e = NaiveEnergyFromDensity(molecs);
        xc_e = XCEnergyFromDensity(molecs,order);

        gauss_e -= 2.0*H2O_gauss_e0;
        naive_e -= 2.0*H2O_e0_naive;
        xc_e -= 2.0.*H2O_e0_XC;

        for i in 1:4
            for ik in 2:(order+1)
                ix = (order+1)*(i-1) + ik;
                if (i < 3)
                    aux_mult_i = ((-1)^(ik-1)) / gamma(1.0 + 2*(ik-1));
                else
                    aux_mult_i = ((-1)^(ik-1)) / gamma(2.0 + 2*(ik-1));
                end

                auxY[i] += aux_mult_i*(gauss_e-naive_e)*xc_e[ix];
                
                for j in 1:4
                    for jk in 2:(order+1)
                        jx = (order+1)*(j-1) + jk;

                        aux_mult_j = 0;
                        if (j < 3)
                            aux_mult_j = ((-1)^(jk-1)) / gamma(1.0 + 2*(jk-1));
                        else
                            aux_mult_j = ((-1)^(jk-1)) / gamma(2.0 + 2*(jk-1));
                        end

                        xc_e_ij = xc_e[ix]*xc_e[jx];
                        aux_mult_ij = aux_mult_i*aux_mult_j;
                        auxM[i,j] += xc_e_ij*aux_mult_ij;
                    end
                end
            end
        end
    end

    return hcat(auxM,auxY);
end

@everywhere function XCCoeffsFromParams(params::Vector)
    xc_coeffs = zeros(typeof(params[1]),4*(order+1));

    for i in 1:2
        i0 = (order+1)*(i-1) + 2;
        i1 = (order+1)*i;

        s = ((-1).^(1:order)) ./ gamma.(1.0 .+ 2.0.*(1:order));
        xc_coeffs[i0:i1] += params[i] .* s;
    end

    for i in 3:4
        i0 = (order+1)*(i-1) + 2;
        i1 = (order+1)*i;

        s = ((-1).^(1:order)) ./ gamma.(2.0 .+ 2.0.*(1:order));
        xc_coeffs[i0:i1] += params[i] .* s;
    end

    return xc_coeffs;
end

aux_mat = pmap(Foo,1:nworkers());

auxM = zeros(Float64,4,4);
auxY = zeros(Float64,4,1);

for i_mat in aux_mat
    global auxM += i_mat[:,1:(end-1)];
    global auxY += i_mat[:,end];
end

fileID = open("ECP_XCCoeffs.txt","w");
for coeff in XCCoeffsFromParams((auxM \ auxY)[:])
    @printf fileID "%22.10E \n" coeff;
end
close(fileID);
