using Profile, BenchmarkTools, Printf
using Optim, Plots, SpecialFunctions, Evolutionary

function trapzz(x::Vector{Float64},y::Vector{Float64})
    coeff = ((x[2]-x[1]))*0.5;

    y = (x.^2.0).*y;
    yy = zeros(Float64,length(x));
    for i in eachindex(yy)
        yy[i] = coeff*(y[1]+y[i]+2.0*sum(y[2:(i-1)]));
    end

    return yy.*(4.0*π);
end

function AtomicNumberElement(atomic_number::Integer)
    # Returns a string of the respective atomic number.
    if atomic_number == 1
        return "H";
    elseif atomic_number == 2
        return "He";
    elseif atomic_number == 3
        return "Li";
    elseif atomic_number == 4
        return "Be";
    elseif atomic_number == 5
        return "B";
    elseif atomic_number == 6
        return "C";
    elseif atomic_number == 7
        return "N";
    elseif atomic_number == 8
        return "O";
    elseif atomic_number == 9
        return "F";
    elseif atomic_number == 10
        return "Ne";
    elseif atomic_number == 11
        return "Na";
    elseif atomic_number == 12
        return "Mg";
    elseif atomic_number == 13
        return "Al";
    elseif atomic_number == 14
        return "Si";
    elseif atomic_number == 15
        return "P";
    elseif atomic_number == 16
        return "S";
    elseif atomic_number == 17
        return "Cl";
    elseif atomic_number == 18
        return "Ar";
    end

    return "X";
end

function GetAtomCharge(file_name::String)
    fileID = open(file_name,"r");

    for i in 1:6
        readline(fileID);
    end

    line = readline(fileID);
    line_splitted = split(line);

    charge = zeros(Float64,1,2);
    charge[1] = parse(Float64,line_splitted[1]);
    charge[2] = parse(Float64,line_splitted[2]);
    return charge;
end

function ReadCubeFile(file_name::String)
    fileID = open(file_name,"r");

    readline(fileID);
    readline(fileID);
    
    grid_data = zeros(Float64,4,4);
    for i in 1:4
        line = readline(fileID);
        line_splitted = split(line);
        for j in 1:4
            grid_data[i,j] = parse(Float64,line_splitted[j]);
        end
    end

    ni = Int(grid_data[2,1]);
    nj = Int(grid_data[3,1]);
    nk = Int(grid_data[4,1]);

    readline(fileID);
    density = zeros(Float64,ni,nj,nk);

    i = 1;
    for line in readlines(fileID)
        line_splitted = split(line);
        for raw_val in line_splitted
            density[i] = parse(Float64,raw_val);
            i = i + 1;
        end
    end
    close(fileID);

    pos_x = zeros(Float64,ni,nj,nk);
    pos_y = zeros(Float64,ni,nj,nk);
    pos_z = zeros(Float64,ni,nj,nk);
    for i in 1:ni
        for j in 1:nj
            for k in 1:nk
                aux_x = grid_data[1,2];
                aux_x += grid_data[2,2]*(i-1);
                aux_x += grid_data[2,3]*(j-1);
                aux_x += grid_data[2,4]*(k-1);

                aux_y = grid_data[1,3];
                aux_y += grid_data[3,2]*(i-1);
                aux_y += grid_data[3,3]*(j-1);
                aux_y += grid_data[3,4]*(k-1);

                aux_z = grid_data[1,4];
                aux_z += grid_data[4,2]*(i-1);
                aux_z += grid_data[4,3]*(j-1);
                aux_z += grid_data[4,4]*(k-1);

                pos_x[i,j,k] = aux_x;
                pos_y[i,j,k] = aux_y;
                pos_z[i,j,k] = aux_z;
            end
        end
    end

    maxl = 10.0;
    num_bins = 1000;

    dl = maxl / num_bins;

    bin_sum = zeros(Float64,num_bins,1);
    sqrt_sum = zeros(Float64,num_bins,1);
    bin_counts = zeros(Float64,num_bins,1);

    for i in eachindex(density)
        px = pos_x[i]
        py = pos_y[i]
        pz = pos_z[i]
        pd = density[i];

        dist = sqrt((px^2.0)+(py^2.0)+(pz^2.0));
        aux_ind = Int32(floor(dist/dl));
        
        if aux_ind < num_bins
            bin_sum[aux_ind+1] += pd;
            sqrt_sum[aux_ind+1] += sqrt(pd);

            bin_counts[aux_ind+1] += 1;
        end
    end

    bin_sum[bin_counts .> 0] ./= bin_counts[bin_counts .> 0];
    sqrt_sum[bin_counts .> 0] ./= bin_counts[bin_counts .> 0];

    data = zeros(Float64,num_bins,5);
    data[:,1] = (0:dl:(maxl-dl));
    data[:,2] = bin_sum;
    data[:,3] = trapzz(data[:,1],data[:,2]);
    data[:,4] = sqrt_sum;
    data[:,5] = trapzz(data[:,1],data[:,4]);

    return data;
end

# vars = 0;
# atomic_number = 3;
# cub_data = zeros(Float64,0,5);

function Foo1(vars::Vector{Float64})
    # Auxiliary function for minimizing the λ value for the electronic 
    # density of an aribitrary atom.
    num_pts = size(cub_data)[1];
    sum1 = 0;
    for i in 1:num_pts
        ρ1i = cub_data[i,2];
        ρ2i = cub_data[i,3];
        di = cub_data[i,1];

        aux_sum1, aux_sum2 = 0,0;
        for j in 1:3
            λ = abs(vars[j]);
            c = vars[3+j];

            aux_sum1 += c*((λ/π)^1.5)*exp(-λ*(di^2.0));

            aux_sum2 += c*erf(di*sqrt(λ));
            aux_sum2 -= 2.0*c*exp(-λ*(di^2.0))*di*sqrt(λ/π);
        end

        sum1 += (aux_sum1-ρ1i)^2.0;
        sum1 += (aux_sum2-ρ2i)^2.0;
    end

    sum2 = 0;
    for i in 1:3
        sum2 += vars[3+i];
    end

    return sum1/num_pts + ((sum2-atomic_number)^2.0);
end

function Foo2(vars::Vector{Float64})
    # Auxiliary function for minimizing the λ value for the electronic 
    # density of an aribitrary atom.
    num_pts = size(cub_data)[1];
    sum1 = 0;
    for i in 1:num_pts
        ρi = cub_data[i,2];
        di = cub_data[i,1];

        aux_sum = 0;

        # core electrons
        for j in 1:6
            λ = abs(vars[j]);
            c = vars[6+j];

            aux_sum += c*((λ/π)^1.5)*exp(-λ*(di^2.0));
        end

        # valence electrons
        for j in 1:3
            λ = abs(ecp_val[j]);
            c = ecp_val[3+j];

            aux_sum += c*((λ/π)^1.5)*exp(-λ*(di^2.0));
        end

        sum1 += (di^(2*pow_val))*((ρi-aux_sum)^2);
    end

    sum2 = 0;

    # core electrons
    for i in 1:6
        sum2 += vars[6+i];
    end

    # valence electrons
    for i in 1:3
        sum2 += ecp_val[3+i];
    end

    return sum1/num_pts + ((sum2-atomic_number)^2.0);
end

function Goo2!(grad::Vector{Float64},vars::Vector{Float64})
    # Gradient of Foo2.
    num_pts = size(cub_data)[1];

    for ii in 1:6
        sum1 = 0;
        grad[6+ii] = 0;

        for i in 1:num_pts
            ρi = cub_data[i,2];
            di = cub_data[i,1];

            aux_sum = 0;

            # core electrons
            for j in 1:6
                λ = abs(vars[j]);
                c = vars[6+j];
    
                aux_sum += c*((λ/π)^1.5)*exp(-λ*(di^2.0));
            end

            # valence electrons
            for j in 1:3
                λ = abs(ecp_val[j]);
                c = ecp_val[3+j];

                aux_sum += c*((λ/π)^1.5)*exp(-λ*(di^2.0));
            end

            λ = abs(vars[ii]);

            aux_mult = ((λ/π)^1.5)*exp(-λ*(di^2.0));
            sum1 -= 2.0*(di^(2*pow_val))*(ρi-aux_sum)*aux_mult;
        end

        sum2 = 0;

        # core electrons
        for i in 1:6
            sum2 += vars[6+i];
        end

        # valence electrons
        for i in 1:3
            sum2 += ecp_val[3+i];
        end

        grad[6+ii] += sum1/num_pts + 2.0*(sum2-atomic_number);
    end
    

    for ii in 1:6
        sum1 = 0;
        grad[ii] = 0;

        for i in 1:num_pts
            ρi = cub_data[i,2];
            di = (cub_data[i,1])^2.0;

            aux_sum = 0;

            # core electrons
            for j in 1:6
                λ = abs(vars[j]);
                c = vars[6+j];
    
                aux_sum += c*((λ/π)^1.5)*exp(-λ*di);
            end

            # valence electrons
            for j in 1:3
                λ = abs(ecp_val[j]);
                c = ecp_val[3+j];
    
                aux_sum += c*((λ/π)^1.5)*exp(-λ*di);
            end

            λ = abs(vars[ii]);
            c = vars[6+ii];

            aux_mult = c*((λ/π)^1.5)*exp(-λ*di)*((1.5/λ)-di);
            sum1 -= 2.0*(di^pow_val)*(ρi-aux_sum)*aux_mult;
        end

        grad[ii] += sum1/num_pts;
    end
end

function GetECPCoeffs()
    # Returns all the coeffients in a convenient matrix form
    coeffs = zeros(Float64,0,6);
    fileID = open("AtomsCoeffs_ECP.txt","r");

    readline(fileID);

    aux_i = 0;
    atom_z = 0;
    for line in readlines(fileID)
        line_splitted = split(line);

        if length(line_splitted) > 6
            aux_i = 1;
            atom_z = parse(Int,line_splitted[2]);
            coeffs = [coeffs; zeros(Float64,1,6)];

            for i in 1:6
                coeffs[end,6*(aux_i-1)+i] = parse(Float64,line_splitted[i+3]);
            end
        else
            aux_i += 1;

            for i in 1:6
                coeffs[end,6*(aux_i-1)+i] = parse(Float64,line_splitted[i]);
            end
        end
    end

    close(fileID);
    return coeffs;
end

ECPCoeffs = GetECPCoeffs();

# rm("AtomsData.txt");
fileID = open("AtomsData.txt","w");
@printf fileID "%10s " "Element";
@printf fileID "%6s " "Z";
@printf fileID "%6s " "Zeff";

for i in 1:3
    @printf fileID "%22s " ("c"*string(i));
    @printf fileID "%22s " ("decay"*string(i));
end
@printf fileID "\n";

ecp_val = zeros(Float64,6);

pow_val = 4.0;
for i in 2:2
    @printf fileID "%10s " AtomicNumberElement(i);
    cub_x = ReadCubeFile("data_base_FullE/"*AtomicNumberElement(i)*"_X.cub");
    cub_y = ReadCubeFile("data_base_FullE/"*AtomicNumberElement(i)*"_Y.cub");
    cub_z = ReadCubeFile("data_base_FullE/"*AtomicNumberElement(i)*"_Z.cub");
    global cub_data = (cub_x + cub_y + cub_z)./3.0;

    charges = GetAtomCharge("data_base_FullE/"*AtomicNumberElement(i)*"_X.cub");
    @printf fileID "%6d " charges[1];
    @printf fileID "%6d " charges[2];

    global atomic_number = charges[2];

    global ecp_val[1:3] = ECPCoeffs[i,2:2:end];
    global ecp_val[4:6] = ECPCoeffs[i,1:2:end];

    init_val = 1.0.*rand(Float64,12);

    opt_val = optimize(Foo2,Goo2!,init_val,ConjugateGradient(),
        Optim.Options(show_trace=true,iterations=6000));

    init_val = Optim.minimizer(opt_val);

    for j in 1:3
        @printf fileID "%22.10f " (Optim.minimizer(opt_val))[6+j];
        @printf fileID "%22.10f " (abs.(Optim.minimizer(opt_val)))[j];
    end
    @printf fileID "\n"

    @printf fileID "%24s " " ";
    for j in 4:6
        @printf fileID "%22.10f " (Optim.minimizer(opt_val))[6+j];
        @printf fileID "%22.10f " (abs.(Optim.minimizer(opt_val)))[j];
    end
    @printf fileID "\n"

    @printf fileID "%24s " " ";
    for j in 1:3
        @printf fileID "%22.10f " ecp_val[3+j];
        @printf fileID "%22.10f " ecp_val[j];
    end
    @printf fileID "\n"

    global vars = Optim.minimizer(opt_val);
    global vars[1:2] = abs.(vars[1:2]);

    println("Done with: "*AtomicNumberElement(i)*".");
end

close(fileID);

plt_1 = 0.0 .* cub_data[:,1];

# core electrons
for i in 1:6
    λ = abs(vars[i]);
    cc = vars[6+i];
    global plt_1 += cc*((λ/π)^1.5)*exp.(-λ.*(cub_data[:,1].^2.0));
end

# valence electrons
for i in 1:3
    λ = abs(ecp_val[i]);
    cc = ecp_val[3+i];
    global plt_1 += cc*((λ/π)^1.5)*exp.(-λ.*(cub_data[:,1].^2.0));
end

p1 = plot(cub_data[:,1],cub_data[:,3]);
plot!(cub_data[:,1],trapzz(cub_data[:,1],plt_1));

p2 = plot(cub_data[:,1],cub_data[:,2]);
plot!(cub_data[:,1],plt_1,ylims=(-0.01,0.05));

plot(p1,p2,layout=(2,1))
