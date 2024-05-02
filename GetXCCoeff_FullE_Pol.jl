using Optim, Distributed, BenchmarkTools

if nworkers() != 11
    addprocs(12-nprocs());
end

@everywhere using ForwardDiff, Random
@everywhere include("GetPotentialFromAngles.jl")

@everywhere all_H2O_H2O_gauss_e = zeros(Float64,10000,1);
@everywhere all_H2O_H2O_charges = zeros(Float64,10000,6);
@everywhere all_H2O_H2O_angs_and_disps = zeros(Float64,10000,6);

@everywhere function ReadGaussianFile(file_name::String, method::String)
    # Reads the energy of the Gaussian file
    energy = 1E30;

    for line in eachline(file_name)
        if method == "DFT"
            if contains(line,"E(RB3LYP)")
                energy = parse(Float64,split(line)[5]);
            end
        elseif method == "CCSD(T)"
            if contains(line,"CCSD(T)= ")
                aux_str = split(line)[end];
                energy = parse(Float64,replace(aux_str, "D" => "E", count = 1));
            end
        end
    end

    return energy;
end

@everywhere function ReadGaussianMullikenCharges(file_name::String)
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

# OH Coord
@everywhere all_angs_and_disps_data = readlines("Training Data/Gaussian Data/FullE/OH Coord/RotationsAndDisplacements.txt");
@everywhere for i in 1:5000
    gauss_file = "Training Data/Gaussian Data/FullE/OH Coord/H2O_H2O_"
    gauss_file *= string(i-1)*".log";

    all_H2O_H2O_gauss_e[i] = ReadGaussianFile(gauss_file,"CCSD(T)");

    local mulliken_charges = ReadGaussianMullikenCharges(gauss_file);
    all_H2O_H2O_charges[i,:] = mulliken_charges;

    angs_and_disps_data = split(all_angs_and_disps_data[i+1]);
    for j in 1:6
        all_H2O_H2O_angs_and_disps[i,j] = parse(Float64,angs_and_disps_data[j]);
    end
end

# OO Coord
@everywhere all_angs_and_disps_data = readlines("Training Data/Gaussian Data/FullE/OO Coord/RotationsAndDisplacements.txt");
@everywhere for i in 1:5000
    gauss_file = "Training Data/Gaussian Data/FullE/OO Coord/H2O_H2O_"
    gauss_file *= string(i-1)*".log";

    all_H2O_H2O_gauss_e[i+5000] = ReadGaussianFile(gauss_file,"CCSD(T)");

    local mulliken_charges = ReadGaussianMullikenCharges(gauss_file);
    all_H2O_H2O_charges[i+5000,:] = mulliken_charges;

    angs_and_disps_data = split(all_angs_and_disps_data[i+1]);
    for j in 1:6
        all_H2O_H2O_angs_and_disps[i+5000,j] = parse(Float64,angs_and_disps_data[j]);
    end
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

@everywhere rand_inds = randperm(10000)[1:20];
@everywhere function RandFooMT(params::Vector)
    # Target function to reach to the get the ideal XC coefficients from the 
    # CCSD(T) water dimer calculations.
    aux_type = typeof(params[1]);
    ret = aux_type(0.0);
    mol_a = ReadMolecule("H2O_fullE_fitted_data.txt",aux_type);

    num_atoms = size(mol_a.atoms_data)[1];
    num_clouds = size(mol_a.cloud_data)[1];
    atoms_per_cloud = ceil(Int,num_clouds/num_atoms);

    xc_coeffs = XCCoeffsFromParams(params);
    for i in rand_inds
        angs_and_disps = all_H2O_H2O_angs_and_disps[i,:];
        
        mol_b = copy(mol_a);
        MoveAndRotateMolec!(mol_b,angs_and_disps);

        molecs = [mol_a,mol_b];
        PolarizeMolecules!(molecs,xc_coeffs);

        model_ρa = 0;
        model_ρb = 0;

        for k in 1:atoms_per_cloud
            new_ρa = mol_a.cloud_data[k:atoms_per_cloud:end,4];
            new_ρa = new_ρa .* mol_a.cloud_data[k:atoms_per_cloud:end,6];
            model_ρa = model_ρa .+ new_ρa;

            new_ρb = mol_b.cloud_data[k:atoms_per_cloud:end,4];
            new_ρb = new_ρb .* mol_b.cloud_data[k:atoms_per_cloud:end,6];
            model_ρb = model_ρb .+ new_ρb;
        end

        model_ρ = vcat(model_ρa,model_ρb);
        target_ρ = vcat(mol_a.atoms_data[:,4],mol_b.atoms_data[:,4]);
        target_ρ -= Vector{aux_type}(all_H2O_H2O_charges[i,:]);
        
        ret += dot(model_ρ - target_ρ,model_ρ - target_ρ);
    end

    return ret / (nworkers()*length(rand_inds));
end

@everywhere function RandFoo(params::Vector)
    return sum(pmap(thread_id -> RandFooMT(params),1:nworkers()));
end

@everywhere function FooMT(params::Vector,thread_id::Int)
    # Target function to reach to the get the ideal XC coefficients from the 
    # CCSD(T) water dimer calculations.
    aux_type = typeof(params[1]);
    ret = aux_type(0.0);
    mol_a = ReadMolecule("H2O_fullE_fitted_data.txt",aux_type);

    num_atoms = size(mol_a.atoms_data)[1];
    num_clouds = size(mol_a.cloud_data)[1];
    atoms_per_cloud = ceil(Int,num_clouds/num_atoms);

    xc_coeffs = XCCoeffsFromParams(params);
    for i in thread_id:nworkers():10000
        angs_and_disps = all_H2O_H2O_angs_and_disps[i,:];
        
        mol_b = copy(mol_a);
        MoveAndRotateMolec!(mol_b,angs_and_disps);

        molecs = [mol_a,mol_b];
        PolarizeMolecules!(molecs,xc_coeffs);

        model_ρa = 0;
        model_ρb = 0;

        for k in 1:atoms_per_cloud
            new_ρa = mol_a.cloud_data[k:atoms_per_cloud:end,4];
            new_ρa = new_ρa .* mol_a.cloud_data[k:atoms_per_cloud:end,6];
            model_ρa = model_ρa .+ new_ρa;

            new_ρb = mol_b.cloud_data[k:atoms_per_cloud:end,4];
            new_ρb = new_ρb .* mol_b.cloud_data[k:atoms_per_cloud:end,6];
            model_ρb = model_ρb .+ new_ρb;
        end

        model_ρ = vcat(model_ρa,model_ρb);
        target_ρ = vcat(mol_a.atoms_data[:,4],mol_b.atoms_data[:,4]);
        target_ρ -= Vector{aux_type}(all_H2O_H2O_charges[i,:]);
        
        ret += dot(model_ρ - target_ρ,model_ρ - target_ρ);
    end

    return ret/10000.0;
end

@everywhere function Foo(params::Vector)
    return sum(pmap(thread_id -> FooMT(params,thread_id),1:nworkers()));
end

@everywhere order = 10;
X0 = zeros(Float64,4);

FullE_result = optimize(RandFoo, X0, NelderMead(), 
    Optim.Options(show_trace = true));
X0 = Optim.minimizer(FullE_result);

FullE_result = optimize(RandFoo,X0, autodiff = :forward,
    LBFGS(), Optim.Options(show_trace = true));
X0 = Optim.minimizer(FullE_result);

FullE_result = optimize(Foo, X0, NelderMead(), 
    Optim.Options(show_trace = true));
X0 = Optim.minimizer(FullE_result);

FullE_result = optimize(Foo,X0, autodiff = :forward,
    LBFGS(), Optim.Options(show_trace = true));
X0 = Optim.minimizer(FullE_result);

xc_coeffs = XCCoeffsFromParams(X0);
fileID = open("FullE_XCCoeffs_Pol.txt","w");
for coeff in xc_coeffs
    @printf fileID "%22.10E \n" coeff;
end
close(fileID);
