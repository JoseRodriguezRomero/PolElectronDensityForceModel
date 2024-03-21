using SpecialFunctions, LinearAlgebra

mutable struct Molecule
    charge::Number;
    atoms_data::Matrix;
    cloud_data::Matrix;
end

Base.copy(mol::Molecule) = Molecule(copy(mol.charge),copy(mol.atoms_data),
    copy(mol.cloud_data));

function Molecule(type = Float64)
    charge = type(0.0);
    atoms_data = zeros(type,0,4);
    cloud_data = zeros(type,0,6);
    return Molecule(charge,atoms_data,cloud_data);
end

function ReadMolecule(molec_data::String, type = Float64)
    fileID = open(molec_data,"r");
    molecule = Molecule(type);

    readline(fileID);

    while true
        line = readline(fileID);
        line_splitted = split(line);

        if length(line_splitted) < 4
            break;
        end

        x = type(parse(Float64,line_splitted[1]));
        y = type(parse(Float64,line_splitted[2]));
        z = type(parse(Float64,line_splitted[3]));
        q = type(parse(Float64,line_splitted[4]));
        
        new_atom = zeros(type,1,4)
        new_atom[1] = x;
        new_atom[2] = y;
        new_atom[3] = z;
        new_atom[4] = q;

        molecule.atoms_data = [molecule.atoms_data; new_atom];
    end

    readline(fileID);
    while true
        line = readline(fileID);
        line_splitted = split(line);

        if length(line_splitted) < 5
            break;
        end

        x = type(parse(Float64,line_splitted[1]));
        y = type(parse(Float64,line_splitted[2]));
        z = type(parse(Float64,line_splitted[3]));
        c = type(parse(Float64,line_splitted[4]));
        λ = type(parse(Float64,line_splitted[5]));

        new_cloud = zeros(type,1,6)
        new_cloud[1] = x;
        new_cloud[2] = y;
        new_cloud[3] = z;
        new_cloud[4] = c;
        new_cloud[5] = λ;
        new_cloud[6] = type(1.0);

        molecule.cloud_data = [molecule.cloud_data; new_cloud];
    end

    close(fileID);
    return molecule;
end

function CenterAtAtomIndex!(molecule::Molecule,index::Int)
    r0 = molecule.atoms_data[index,1:3];
    molecule.atoms_data[:,1:3] .-= r0';
    molecule.cloud_data[:,1:3] .-= r0';
end

function MoveAndRotateMolec!(molecule::Molecule,angs_and_disp::Vector)
    θx = angs_and_disp[1];
    θy = angs_and_disp[2];
    θz = angs_and_disp[3];
    dr = angs_and_disp[4:6];

    aux_type = typeof(angs_and_disp[1,1]);

    rot_x = zeros(aux_type,3,3);
    rot_x[1,1] = 1.0;
    rot_x[2,2] = cos(θx);
    rot_x[2,3] = -sin(θx);
    rot_x[3,2] = sin(θx);
    rot_x[3,3] = cos(θx);

    rot_y = zeros(aux_type,3,3);
    rot_y[1,1] = cos(θy);
    rot_y[1,3] = sin(θy);
    rot_y[2,2] = 1.0;
    rot_y[3,1] = -sin(θy);
    rot_y[3,3] = cos(θy);

    rot_z = zeros(aux_type,3,3);
    rot_z[1,1] = cos(θz);
    rot_z[1,2] = -sin(θz);
    rot_z[2,1] = sin(θz);
    rot_z[2,2] = cos(θz);
    rot_z[3,3] = 1.0;

    rot = (rot_x*rot_y*rot_z)';

    molecule.atoms_data[:,1:3] = molecule.atoms_data[:,1:3]*rot;
    molecule.atoms_data[:,1:3] .+= dr';

    molecule.cloud_data[:,1:3] = molecule.cloud_data[:,1:3]*rot;
    molecule.cloud_data[:,1:3] .+= dr';
end

function XCOrder0(λ::Vector,l::Vector)
    if abs(l[1]) < 1.0E-6
        return 2*sqrt(λ[1]/π);
    else
        return erf(l[1]sqrt(λ[1]))/l[1];
    end
end

function XCOrder1(λ::Vector)
    return -((4*λ[1]^(3/2))/sqrt(π))
end

function XCOrder2(λ::Vector,l::Vector)
    return -((sqrt(λ[1])*(16*l[2]*λ[3]-24*λ[2]))/sqrt(π));
end

function XCOrder3(λ::Vector,l::Vector)
    return -((sqrt(λ[1])*(64*l[4]*λ[5]-320*l[2]*λ[4]+240*λ[3]))/sqrt(π));
end

function XCOrder4(λ::Vector,l::Vector)
    return -((sqrt(λ[1])*(256*l[6]*λ[7]-2688*l[4]*λ[6]+6720*l[2]*λ[5]-
        3360*λ[4]))/sqrt(π));
end

function XCOrder5(λ::Vector,l::Vector)
    return -((sqrt(λ[1])*(1024*l[8]*λ[9]-18432*l[6]*λ[8]+96768*l[4]*λ[7]-
        161280*l[2]*λ[6]+60480*λ[5]))/sqrt(π));
end

function XCOrder6(λ::Vector,l::Vector)
    return -((sqrt(λ[1])*(4096*l[10]*λ[11]-112640*l[8]*λ[10]+1013760*l[6]*λ[9]-
        3548160*l[4]*λ[8]+4435200*l[2]*λ[7]-1330560*λ[6]))/sqrt(π));
end

function XCOrder7(λ::Vector,l::Vector)
    return -((sqrt(λ[1])*(16384*l[12]*λ[13]-638976*l[10]*λ[12]+8785920*l[8]*
        λ[11]-52715520*l[6]*λ[10]+138378240*l[4]*λ[9]-138378240*l[2]*λ[8]+
        34594560*λ[7]))/sqrt(π));
end

function XCOrder8(λ::Vector,l::Vector)
    return -((sqrt(λ[1])*(65536*l[14]*λ[15]-3440640*l[12]*λ[14]+67092480*l[10]*
        λ[13]-615014400*l[8]*λ[12]+2767564800*l[6]*λ[11]-5811886080*l[4]*λ[10]+
        4843238400*l[2]*λ[9]-1037836800*λ[8]))/sqrt(π));
end

function XCOrder9(λ::Vector,l::Vector)
    return -((sqrt(λ[1])*(262144*l[16]*λ[17]-17825792*l[14]*λ[16]+467927040*
        l[12]*λ[15]-6083051520*l[10]*λ[14]+41820979200*l[8]*λ[13]-150555525120*
        l[6]*λ[12]+263472168960*l[4]*λ[11]-188194406400*l[2]*λ[10]+35286451200*
        λ[9]))/sqrt(π));
end

function XCOrder10(λ::Vector,l::Vector)
    return -((sqrt(λ[1])*(1048576*l[18]*λ[19]-89653248*l[16]*λ[18]+3048210432*
        l[14]*λ[17]-53343682560*l[12]*λ[16]+520100904960*l[10]*λ[15]-
        2860554977280*l[8]*λ[14]+8581664931840*l[6]*λ[13]-12872497397760*l[4]*
        λ[12]+8045310873600*l[2]*λ[11]-1340885145600*λ[10]))/sqrt(π));
end

function XCOrder11(λ::Vector,l::Vector)
    return -((sqrt(λ[1])*(4194304*l[20]*λ[21]-440401920*l[18]*λ[20]+18827182080*
        l[16]*λ[19]-426749460480*l[14]*λ[18]+5601086668800*l[12]*λ[17]-
        43688476016640*l[10]*λ[16]+200238848409600*l[8]*λ[15]-514899895910400*
        l[6]*λ[14]+675806113382400*l[4]*λ[13]-375447840768000*l[2]*
        λ[12]+56317176115200*λ[11]))/sqrt(π));
end

function XCOrderD00(λ::Vector,l::Vector)
    if abs(l[1]) < 1.0E-8
        return 0.0;
    else
        return 2*l[1]sqrt(λ[1])/(sqrt(π)*l[2]);
    end
end

function XCOrderD01(λ::Vector,l::Vector)
    if abs(l[1]) < 1.0E-8
        return 0.0;
    else
        return -sqrt(π)*erf(l[1]sqrt(λ[1]))/(sqrt(π)*l[2]);;
    end
end

function XCOrderD1(λ::Vector,l::Vector)
    return (8*l[1]*λ[1]^(5/2))/sqrt(π);
end

function XCOrderD2(λ::Vector,l::Vector)
    return (sqrt(λ[1])*(32*l[3]*λ[4]-80*l[1]*λ[3]))/sqrt(π);
end

function XCOrderD3(λ::Vector,l::Vector)
    return (sqrt(λ[1])*(128*l[5]*λ[6]-896*l[3]*λ[5]+1120*l[1]*λ[4]))/sqrt(π);
end

function XCOrderD4(λ::Vector,l::Vector)
    return (sqrt(λ[1])*(512*l[7]*λ[8]-6912*l[5]*λ[7]+24192*l[3]*λ[6]-20160*l[1]*
        λ[5]))/sqrt(π);
end

function XCOrderD5(λ::Vector,l::Vector)
    return (sqrt(λ[1])*(2048*l[9]*λ[10]-45056*l[7]*λ[9]+304128*l[5]*λ[8]-709632*
        l[3]*λ[7]+443520*l[1]*λ[6]))/sqrt(π);
end

function XCOrderD6(λ::Vector,l::Vector)
    return (sqrt(λ[1])*(8192*l[11]*λ[12]-266240*l[9]*λ[11]+2928640*l[7]*λ[10]-
        13178880*l[5]*λ[9]+23063040*l[3]*λ[8]-11531520*l[1]*λ[7]))/sqrt(π);
end

function XCOrderD7(λ::Vector,l::Vector)
    return (sqrt(λ[1])*(32768*l[13]*λ[14]-1474560*l[11]*λ[13]+23961600*l[9]*
        λ[12]-175718400*l[7]*λ[11]+593049600*l[5]*λ[10]-830269440*l[3]*λ[9]+
        345945600*l[1]*λ[8]))/sqrt(π);
end

function XCOrderD8(λ::Vector,l::Vector)
    return (sqrt(λ[1])*(131072*l[15]*λ[16]-7798784*l[13]*λ[15]+175472640*l[11]*
        λ[14]-1900953600*l[9]*λ[13]+10455244800*l[7]*λ[12]-28229160960*l[5]*
        λ[11]+32934021120*l[3]*λ[10]-11762150400*l[1]*λ[9]))/sqrt(π);
end

function XCOrderD9(λ::Vector,l::Vector)
    return (sqrt(λ[1])*(524288*l[17]*λ[18]-39845888*l[15]*λ[17]+1185415168*
        l[13]*λ[16]-17781227520*l[11]*λ[15]+144472473600*l[9]*λ[14]-
        635678883840*l[7]*λ[13]+1430277488640*l[5]*λ[12]-1430277488640*l[3]*
        λ[11]+446961715200*l[1]*λ[10]))/sqrt(π);
end

function XCOrderD10(λ::Vector,l::Vector)
    return (sqrt(λ[1])*(2097152*l[19]*λ[20]-198180864*l[17]*λ[19]+7530872832*
        l[15]*λ[18]-149362311168*l[13]*λ[17]+1680326000640*l[11]*λ[16]-
        10922119004160*l[9]*λ[15]+40047769681920*l[7]*λ[14]-77234984386560*l[5]*
        λ[13]+67580611338240*l[3]*λ[12]-18772392038400*l[1]*λ[11]))/sqrt(π);
end

function XCOrderD11(λ::Vector,l::Vector)
    return (sqrt(λ[1])*(8388608*l[21]*λ[22]-964689920*l[19]*λ[21]+45581598720*
        l[17]*λ[20]-1154733834240*l[15]*λ[19]+17176665784320*l[13]*λ[18]-
        154589992058880*l[11]*λ[17]+837362456985600*l[9]*λ[16]-2631710579097600*
        l[7]*λ[15]+4441011602227200*l[5]*λ[14]-3454120135065600*l[3]*λ[13]+
        863530033766400*l[1]*λ[12]))/sqrt(π);
end

function PolarizeMolecules!(molecules::Vector,xc_coeffs::Vector)
    # Polarizes the molecules by changing their polarization coeffients (evert 
    # sixth entry molecule.cloud_data[:,6]) for being all ones (non-polarized) 
    # to something else.
    order = ceil(Int,length(xc_coeffs)/4) - 1;
    num_molecules = length(molecules);

    aux_type = typeof(xc_coeffs[1]);

    tot_num_vars = 0;
    molec_ind_base = [];
    for molecule in molecules
        if isempty(molec_ind_base)
            push!(molec_ind_base,1);
            tot_num_vars += length(molecule.atoms_data[:,1]);
        else
            aux_num_atoms = size(molecule.atoms_data)[1];
            push!(molec_ind_base,molec_ind_base[end]+aux_num_atoms);

            tot_num_vars += aux_num_atoms;
        end
    end

    aux_M = zeros(aux_type,tot_num_vars+1,tot_num_vars+1);
    aux_Y = zeros(aux_type,tot_num_vars+1,1);

    for ii in eachindex(molecules)
        molecule1 = molecules[ii];
        num_atoms1 = size(molecule1.atoms_data)[1];
        num_clouds1 = size(molecule1.cloud_data)[1];
        clouds_per_atom = ceil(Int,num_clouds1/num_atoms1);

        frozen_clouds1 = [];
        valence_clouds1 = [];
        if clouds_per_atom > 3
            for i in 1:6
                frozen_clouds1 = vcat(frozen_clouds1,i:clouds_per_atom:num_clouds1);
            end

            for i in 7:9
                valence_clouds1 = vcat(valence_clouds1,i:clouds_per_atom:num_clouds1);
            end
        else
            frozen_clouds1 = [];
            valence_clouds1 = 1:num_clouds1;
        end

        # Lagrange multiplier to keep the number of electrons constant 
        # in the new fictitious electron density.
        aux_Y[end] -= molecules[ii].charge;
        
        ii0 = molec_ind_base[ii];
        if clouds_per_atom > 3
            for i in 1:num_atoms1
                for j in (clouds_per_atom-2):clouds_per_atom
                    ij = clouds_per_atom*(i-1) + j;
                    aux_M[end,ii0+i-1] += molecule1.cloud_data[ij,4];
                    aux_M[ii0+i-1,end] += molecule1.cloud_data[ij,4];
                    aux_Y[end] += molecule1.cloud_data[ij,4];
                end
            end
        else
            for i in 1:num_atoms1
                for j in 1:clouds_per_atom
                    ij = clouds_per_atom*(i-1) + j;
                    aux_M[end,ii0+i-1] += molecule1.cloud_data[ij,4];
                    aux_M[ii0+i-1,end] += molecule1.cloud_data[ij,4];
                    aux_Y[end] += molecule1.cloud_data[ij,4];
                end
            end
        end

        for jj in eachindex(molecules)
            molecule2 = molecules[jj];
            num_atoms2 = size(molecule2.atoms_data)[1];
            num_clouds2 = size(molecule2.cloud_data)[1];

            frozen_clouds2 = [];
            valence_clouds2 = [];
            if clouds_per_atom > 3
                for i in 1:6
                    frozen_clouds2 = vcat(frozen_clouds2,i:clouds_per_atom:num_clouds2);
                end

                for i in 7:9
                    valence_clouds2 = vcat(valence_clouds2,i:clouds_per_atom:num_clouds2);
                end
            else
                frozen_clouds2 = [];
                valence_clouds2 = 1:num_clouds2;
            end

            # valence_cloud - atom
            for i in valence_clouds1
                for j in 1:num_atoms2
                    c1 = molecule1.cloud_data[i,4];
                    λ1 = molecule1.cloud_data[i,5];
                    q2 = molecule2.atoms_data[j,4];

                    if c1 == 0
                        continue;
                    end

                    aux_dist = molecule1.cloud_data[i,1:3];
                    aux_dist -= molecule2.atoms_data[j,1:3];
                    aux_dist = norm(aux_dist);

                    C1 = copy(c1);
                    c1 *= exp(-λ1*(aux_dist^2.0));

                    aux_dist = zeros(aux_type,22) .+ aux_dist;
                    aux_dist = aux_dist .^ collect(1:22);

                    λ1 = zeros(aux_type,22) .+ λ1;
                    λ1 = λ1 .^ collect(1:22);

                    ii0 = molec_ind_base[ii] + ceil(Int,i/clouds_per_atom) - 1;

                    # Naive contribution
                    aux_Y[ii0] -= C1*q2*XCOrder0(λ1,aux_dist);

                    # XC contributions
                    for k_order in 1:(order+1)
                        xc_coeff_1 = xc_coeffs[0*(order+1)+k_order];
                        xc_coeff_2 = xc_coeffs[2*(order+1)+k_order];
                        if k_order == 1
                            aux_Y[ii0] += C1*q2*xc_coeff_1*XCOrder0(λ1,aux_dist);
                            aux_Y[ii0] += C1*q2*xc_coeff_2*XCOrderD01(λ1,aux_dist);
                            aux_Y[ii0] += c1*q2*xc_coeff_2*XCOrderD00(λ1,aux_dist);
                        elseif k_order == 2
                            aux_Y[ii0] += c1*q2*xc_coeff_1*XCOrder1(λ1);
                            aux_Y[ii0] += c1*q2*xc_coeff_2*XCOrderD1(λ1,aux_dist);
                        elseif k_order == 3
                            aux_Y[ii0] += c1*q2*xc_coeff_1*XCOrder2(λ1,aux_dist);
                            aux_Y[ii0] += c1*q2*xc_coeff_2*XCOrderD2(λ1,aux_dist);
                        elseif k_order == 4
                            aux_Y[ii0] += c1*q2*xc_coeff_1*XCOrder3(λ1,aux_dist);
                            aux_Y[ii0] += c1*q2*xc_coeff_2*XCOrderD3(λ1,aux_dist);
                        elseif k_order == 5
                            aux_Y[ii0] += c1*q2*xc_coeff_1*XCOrder4(λ1,aux_dist);
                            aux_Y[ii0] += c1*q2*xc_coeff_2*XCOrderD4(λ1,aux_dist);
                        elseif k_order == 6
                            aux_Y[ii0] += c1*q2*xc_coeff_1*XCOrder5(λ1,aux_dist);
                            aux_Y[ii0] += c1*q2*xc_coeff_2*XCOrderD5(λ1,aux_dist);
                        elseif k_order == 7
                            aux_Y[ii0] += c1*q2*xc_coeff_1*XCOrder6(λ1,aux_dist);
                            aux_Y[ii0] += c1*q2*xc_coeff_2*XCOrderD6(λ1,aux_dist);
                        elseif k_order == 8
                            aux_Y[ii0] += c1*q2*xc_coeff_1*XCOrder7(λ1,aux_dist);
                            aux_Y[ii0] += c1*q2*xc_coeff_2*XCOrderD7(λ1,aux_dist);
                        elseif k_order == 9
                            aux_Y[ii0] += c1*q2*xc_coeff_1*XCOrder8(λ1,aux_dist);
                            aux_Y[ii0] += c1*q2*xc_coeff_2*XCOrderD8(λ1,aux_dist);
                        elseif k_order == 10
                            aux_Y[ii0] += c1*q2*xc_coeff_1*XCOrder9(λ1,aux_dist);
                            aux_Y[ii0] += c1*q2*xc_coeff_2*XCOrderD9(λ1,aux_dist);
                        elseif k_order == 11
                            aux_Y[ii0] += c1*q2*xc_coeff_1*XCOrder10(λ1,aux_dist);
                            aux_Y[ii0] += c1*q2*xc_coeff_2*XCOrderD10(λ1,aux_dist);
                        end
                    end
                end
            end

            # valence_cloud - frozen_cloud
            for i in valence_clouds1
                for j in frozen_clouds2
                    c1 = molecule1.cloud_data[i,4];
                    λ1 = molecule1.cloud_data[i,5];
                    c2 = molecule2.cloud_data[j,4];
                    λ2 = molecule2.cloud_data[j,5];

                    if c1 == 0
                        continue;
                    end

                    λ = (λ1*λ2)/(λ1+λ2);

                    aux_dist = molecule1.cloud_data[i,1:3];
                    aux_dist -= molecule2.cloud_data[j,1:3];
                    aux_dist = norm(aux_dist);

                    C1 = copy(c1);
                    c1 *= exp(-λ*(aux_dist^2.0));

                    aux_dist = zeros(aux_type,22) .+ aux_dist;
                    aux_dist = aux_dist .^ collect(1:22);

                    λ = zeros(aux_type,22) .+ λ;
                    λ = λ .^ collect(1:22);

                    ii0 = molec_ind_base[ii] + ceil(Int,i/clouds_per_atom) - 1;

                    # Naive contribution
                    aux_Y[ii0] += C1*c2*XCOrder0(λ,aux_dist);

                    # XC contributions
                    for k_order in 1:(order+1)
                        xc_coeff_1 = xc_coeffs[1*(order+1)+k_order];
                        xc_coeff_2 = xc_coeffs[3*(order+1)+k_order];
                        if k_order == 1
                            aux_Y[ii0] += C1*c2*xc_coeff_1*XCOrder0(λ,aux_dist);
                            aux_Y[ii0] += C1*c2*xc_coeff_2*XCOrderD01(λ,aux_dist);
                            aux_Y[ii0] += c1*c2*xc_coeff_2*XCOrderD00(λ,aux_dist);
                        elseif k_order == 2
                            aux_Y[ii0] += c1*c2*xc_coeff_1*XCOrder1(λ);
                            aux_Y[ii0] += c1*c2*xc_coeff_2*XCOrderD1(λ,aux_dist);
                        elseif k_order == 3
                            aux_Y[ii0] += c1*c2*xc_coeff_1*XCOrder2(λ,aux_dist);
                            aux_Y[ii0] += c1*c2*xc_coeff_2*XCOrderD2(λ,aux_dist);
                        elseif k_order == 4
                            aux_Y[ii0] += c1*c2*xc_coeff_1*XCOrder3(λ,aux_dist);
                            aux_Y[ii0] += c1*c2*xc_coeff_2*XCOrderD3(λ,aux_dist);
                        elseif k_order == 5
                            aux_Y[ii0] += c1*c2*xc_coeff_1*XCOrder4(λ,aux_dist);
                            aux_Y[ii0] += c1*c2*xc_coeff_2*XCOrderD4(λ,aux_dist);
                        elseif k_order == 6
                            aux_Y[ii0] += c1*c2*xc_coeff_1*XCOrder5(λ,aux_dist);
                            aux_Y[ii0] += c1*c2*xc_coeff_2*XCOrderD5(λ,aux_dist);
                        elseif k_order == 7
                            aux_Y[ii0] += c1*c2*xc_coeff_1*XCOrder6(λ,aux_dist);
                            aux_Y[ii0] += c1*c2*xc_coeff_2*XCOrderD6(λ,aux_dist);
                        elseif k_order == 8
                            aux_Y[ii0] += c1*c2*xc_coeff_1*XCOrder7(λ,aux_dist);
                            aux_Y[ii0] += c1*c2*xc_coeff_2*XCOrderD7(λ,aux_dist);
                        elseif k_order == 9
                            aux_Y[ii0] += c1*c2*xc_coeff_1*XCOrder8(λ,aux_dist);
                            aux_Y[ii0] += c1*c2*xc_coeff_2*XCOrderD8(λ,aux_dist);
                        elseif k_order == 10
                            aux_Y[ii0] += c1*c2*xc_coeff_1*XCOrder9(λ,aux_dist);
                            aux_Y[ii0] += c1*c2*xc_coeff_2*XCOrderD9(λ,aux_dist);
                        elseif k_order == 11
                            aux_Y[ii0] += c1*c2*xc_coeff_1*XCOrder10(λ,aux_dist);
                            aux_Y[ii0] += c1*c2*xc_coeff_2*XCOrderD10(λ,aux_dist);
                        end
                    end
                end
            end

            # valence_cloud - valence_cloud
            for i in valence_clouds1
                for j in valence_clouds2
                    c1 = molecule1.cloud_data[i,4];
                    λ1 = molecule1.cloud_data[i,5];
                    c2 = molecule2.cloud_data[j,4];
                    λ2 = molecule2.cloud_data[j,5];

                    if (c1 == 0) || (c2 == 0)
                        continue;
                    end

                    λ = (λ1*λ2)/(λ1+λ2);

                    aux_dist = molecule1.cloud_data[i,1:3];
                    aux_dist -= molecule2.cloud_data[j,1:3];
                    aux_dist = norm(aux_dist);

                    ii0 = molec_ind_base[ii] + ceil(Int,i/clouds_per_atom) - 1;
                    jj0 = molec_ind_base[jj] + ceil(Int,j/clouds_per_atom) - 1;

                    if (ii == jj) && (ii0 == jj0)
                        c1 *= 2.0;
                    end

                    C1 = copy(c1);
                    c1 *= exp(-λ*(aux_dist^2.0));

                    aux_dist = zeros(aux_type,22) .+ aux_dist;
                    aux_dist = aux_dist .^ collect(1:22);

                    λ = zeros(aux_type,22) .+ λ;
                    λ = λ .^ collect(1:22);

                    # Naive contribution
                    aux_M[ii0,jj0] += C1*c2*XCOrder0(λ,aux_dist);

                    # XC contributions
                    for k_order in 1:(order+1)
                        xc_coeff_1 = xc_coeffs[1*(order+1)+k_order];
                        xc_coeff_2 = xc_coeffs[3*(order+1)+k_order];
                        if k_order == 1
                            aux_M[ii0,jj0] += C1*c2*xc_coeff_1*XCOrder0(λ,aux_dist);
                            aux_M[ii0,jj0] += C1*c2*xc_coeff_2*XCOrderD01(λ,aux_dist);
                            aux_M[ii0,jj0] += c1*c2*xc_coeff_2*XCOrderD00(λ,aux_dist);
                        elseif k_order == 2
                            aux_M[ii0,jj0] += c1*c2*xc_coeff_1*XCOrder1(λ);
                            aux_M[ii0,jj0] += c1*c2*xc_coeff_2*XCOrderD1(λ,aux_dist);
                        elseif k_order == 3
                            aux_M[ii0,jj0] += c1*c2*xc_coeff_1*XCOrder2(λ,aux_dist);
                            aux_M[ii0,jj0] += c1*c2*xc_coeff_2*XCOrderD2(λ,aux_dist);
                        elseif k_order == 4
                            aux_M[ii0,jj0] += c1*c2*xc_coeff_1*XCOrder3(λ,aux_dist);
                            aux_M[ii0,jj0] += c1*c2*xc_coeff_2*XCOrderD3(λ,aux_dist);
                        elseif k_order == 5
                            aux_M[ii0,jj0] += c1*c2*xc_coeff_1*XCOrder4(λ,aux_dist);
                            aux_M[ii0,jj0] += c1*c2*xc_coeff_2*XCOrderD4(λ,aux_dist);
                        elseif k_order == 6
                            aux_M[ii0,jj0] += c1*c2*xc_coeff_1*XCOrder5(λ,aux_dist);
                            aux_M[ii0,jj0] += c1*c2*xc_coeff_2*XCOrderD5(λ,aux_dist);
                        elseif k_order == 7
                            aux_M[ii0,jj0] += c1*c2*xc_coeff_1*XCOrder6(λ,aux_dist);
                            aux_M[ii0,jj0] += c1*c2*xc_coeff_2*XCOrderD6(λ,aux_dist);
                        elseif k_order == 8
                            aux_M[ii0,jj0] += c1*c2*xc_coeff_1*XCOrder7(λ,aux_dist);
                            aux_M[ii0,jj0] += c1*c2*xc_coeff_2*XCOrderD7(λ,aux_dist);
                        elseif k_order == 9
                            aux_M[ii0,jj0] += c1*c2*xc_coeff_1*XCOrder8(λ,aux_dist);
                            aux_M[ii0,jj0] += c1*c2*xc_coeff_2*XCOrderD8(λ,aux_dist);
                        elseif k_order == 10
                            aux_M[ii0,jj0] += c1*c2*xc_coeff_1*XCOrder9(λ,aux_dist);
                            aux_M[ii0,jj0] += c1*c2*xc_coeff_2*XCOrderD9(λ,aux_dist);
                        elseif k_order == 11
                            aux_M[ii0,jj0] += c1*c2*xc_coeff_1*XCOrder10(λ,aux_dist);
                            aux_M[ii0,jj0] += c1*c2*xc_coeff_2*XCOrderD10(λ,aux_dist);
                        end
                    end
                end
            end
        end
    end

    aux_Y[1:(end-1)] .*= -1.0;
    aux_X = aux_M \ aux_Y;

    # display(hcat(aux_M,aux_Y));

    for ii in eachindex(molecules)
        ii0 = molec_ind_base[ii];
        molecule1 = molecules[ii];
        num_atoms1 = size(molecule1.atoms_data)[1];
        num_clouds1 = size(molecule1.cloud_data)[1];
        clouds_per_atom = ceil(Int,num_clouds1/num_atoms1);

        if clouds_per_atom == 3
            χ = aux_X[ii0:(ii0+num_atoms1-1)];
            for i in 1:clouds_per_atom
                molecules[ii].cloud_data[i:clouds_per_atom:end,6] = χ;
            end
        else
            χ = aux_X[ii0:(ii0+num_atoms1-1)];
            for i in 7:clouds_per_atom
                molecules[ii].cloud_data[i:clouds_per_atom:end,6] = χ;
            end
        end 
    end
end

function NaiveEnergyFromDensity(molecules::Vector{Molecule})
    # Returns the surface energy from the fitted electron clouds and nuclei
    # positions using the polarized naive model.
    energy = 0;
    num_molecules = length(molecules);

    aux_type = typeof(molecules[1].atoms_data[1,1]);

    for ii in 1:num_molecules
        molecule1 = molecules[ii];
        num_atoms1 = size(molecule1.atoms_data)[1];
        num_clouds1 = size(molecule1.cloud_data)[1];

        # Intramolecular interactions
        # nuclei-nuclei
        for i in 1:num_atoms1
            for j in (i+1):num_atoms1
                q1 = molecule1.atoms_data[i,4];
                q2 = molecule1.atoms_data[j,4];

                aux_dist = molecule1.atoms_data[i,1:3];
                aux_dist -= molecule1.atoms_data[j,1:3];
                aux_dist = norm(aux_dist);

                energy += q1*q2/aux_dist;
            end
        end

        # nuclei-cloud
        for i in 1:num_atoms1
            for j in 1:num_clouds1
                q1 = molecule1.atoms_data[i,4];
                c2 = molecule1.cloud_data[j,4];
                λ2 = molecule1.cloud_data[j,5];

                if c2 == 0
                    continue;
                end

                # polarization
                c2 *= molecule1.cloud_data[j,6];

                aux_dist = molecule1.atoms_data[i,1:3];
                aux_dist -= molecule1.cloud_data[j,1:3];
                aux_dist = norm(aux_dist);

                aux_dist = zeros(aux_type,1) .+ aux_dist;
                aux_dist = aux_dist .^collect(1:1);

                λ2 = zeros(aux_type,1) .+ λ2;
                λ2 = λ2 .^ collect(1:1);

                energy -= q1*c2*XCOrder0(λ2,aux_dist);
            end
        end

        # cloud-cloud
        for i in 1:num_clouds1
            for j in i:num_clouds1
                c1 = molecule1.cloud_data[i,4];
                λ1 = molecule1.cloud_data[i,5];
                c2 = molecule1.cloud_data[j,4];
                λ2 = molecule1.cloud_data[j,5];

                if (c1 == 0) || (c2 == 0)
                    continue;
                end

                λ = (λ1*λ2)/(λ1+λ2);

                # polarization
                c1 *= molecule1.cloud_data[i,6];
                c2 *= molecule1.cloud_data[j,6];

                aux_dist = molecule1.cloud_data[i,1:3];
                aux_dist -= molecule1.cloud_data[j,1:3];
                aux_dist = norm(aux_dist);

                aux_dist = zeros(aux_type,1) .+ aux_dist;
                aux_dist = aux_dist .^collect(1:1);

                λ = zeros(aux_type,1) .+ λ;
                λ = λ .^ collect(1:1);

                energy += c1*c2*XCOrder0(λ,aux_dist);
            end
        end

        # Intermolecular interactions
        for jj in (ii+1):num_molecules
            molecule2 = molecules[jj];
            num_atoms2 = size(molecule2.atoms_data)[1];
            num_clouds2 = size(molecule2.cloud_data)[1];

            # nuclei-nuclei
            for i in 1:num_atoms1
                for j in 1:num_atoms2
                    q1 = molecule1.atoms_data[i,4];
                    q2 = molecule2.atoms_data[j,4];
    
                    aux_dist = molecule1.atoms_data[i,1:3];
                    aux_dist -= molecule2.atoms_data[j,1:3];
                    aux_dist = norm(aux_dist);
    
                    energy += q1*q2/aux_dist;
                end
            end

            # nuclei-cloud
            for i in 1:num_atoms1
                for j in 1:num_clouds2
                    q1 = molecule1.atoms_data[i,4];
                    c2 = molecule2.cloud_data[j,4];
                    λ2 = molecule2.cloud_data[j,5];

                    if c2 == 0
                        continue;
                    end

                    # polarization
                    c2 *= molecule2.cloud_data[j,6];

                    aux_dist = molecule1.atoms_data[i,1:3];
                    aux_dist -= molecule2.cloud_data[j,1:3];
                    aux_dist = norm(aux_dist);

                    aux_dist = zeros(aux_type,1) .+ aux_dist;
                    aux_dist = aux_dist .^ collect(1:1);

                    λ2 = zeros(aux_type,1) .+ λ2;
                    λ2 = λ2 .^ collect(1:1);

                    energy -= q1*c2*XCOrder0(λ2,aux_dist);
                end
            end

            for i in 1:num_clouds1
                for j in 1:num_atoms2
                    c1 = molecule1.cloud_data[i,4];
                    λ1 = molecule1.cloud_data[i,5];
                    q2 = molecule2.atoms_data[j,4];

                    if c1 == 0
                        continue;
                    end

                    # polarization
                    c1 *= molecule1.cloud_data[i,6];

                    aux_dist = molecule1.cloud_data[i,1:3];
                    aux_dist -= molecule2.atoms_data[j,1:3];
                    aux_dist = norm(aux_dist);

                    aux_dist = zeros(aux_type,1) .+ aux_dist;
                    aux_dist = aux_dist .^ collect(1:1);

                    λ1 = zeros(aux_type,1) .+ λ1;
                    λ1 = λ1 .^ collect(1:1);

                    energy -= c1*q2*XCOrder0(λ1,aux_dist);
                end
            end

            # cloud-cloud
            for i in 1:num_clouds1
                for j in 1:num_clouds2
                    c1 = molecule1.cloud_data[i,4];
                    λ1 = molecule1.cloud_data[i,5];
                    c2 = molecule2.cloud_data[j,4];
                    λ2 = molecule2.cloud_data[j,5];

                    if (c1 == 0) || (c2 == 0)
                        continue;
                    end

                    λ = (λ1*λ2)/(λ1+λ2);

                    # polarization
                    c1 *= molecule1.cloud_data[i,6];
                    c2 *= molecule2.cloud_data[j,6];

                    aux_dist = molecule1.cloud_data[i,1:3];
                    aux_dist -= molecule2.cloud_data[j,1:3];
                    aux_dist = norm(aux_dist);

                    aux_dist = zeros(aux_type,1) .+ aux_dist;
                    aux_dist = aux_dist .^ collect(1:1);

                    λ = zeros(aux_type,1) .+ λ;
                    λ = λ .^ collect(1:1);

                    energy += c1*c2*XCOrder0(λ,aux_dist);
                end
            end
        end
    end

    return energy;
end

function XCEnergyFromDensity(molecules::Vector{Molecule},order::Int)
    # Returns the surface XC energy from the fitted polarized electron clouds 
    # and nuclei positions.
    aux_type = typeof(molecules[1].atoms_data[1,1]);
    energies = zeros(aux_type,1,4*(order+1));
    num_molecules = length(molecules);

    for ii in 1:num_molecules
        molecule1 = molecules[ii];
        num_atoms1 = size(molecule1.atoms_data)[1];
        num_clouds1 = size(molecule1.cloud_data)[1];

        # Intramolecular interactions
        # nuclei-cloud 
        for i in 1:num_atoms1
            for j in 1:num_clouds1
                q1 = molecule1.atoms_data[i,4];
                c2 = molecule1.cloud_data[j,4];
                λ2 = molecule1.cloud_data[j,5];

                if c2 == 0
                    continue;
                end

                # polarization
                c2 *= molecule1.cloud_data[j,6];

                aux_dist = molecule1.atoms_data[i,1:3];
                aux_dist -= molecule1.cloud_data[j,1:3];
                aux_dist = norm(aux_dist);

                C2 = copy(c2);
                c2 *= exp(-λ2*(aux_dist^2.0));

                aux_dist = zeros(aux_type,22) .+ aux_dist;
                aux_dist = aux_dist .^ collect(1:22);

                λ2 = zeros(aux_type,22) .+ λ2;
                λ2 = λ2 .^ collect(1:22);
                for k_order in 1:(order+1)
                    if k_order == 1
                        energies[0*(order+1)+k_order] += q1*C2*XCOrder0(λ2,aux_dist);
                        energies[2*(order+1)+k_order] += q1*C2*XCOrderD01(λ2,aux_dist);
                        energies[2*(order+1)+k_order] += q1*c2*XCOrderD00(λ2,aux_dist);
                    elseif k_order == 2
                        energies[0*(order+1)+k_order] += q1*c2*XCOrder1(λ2);
                        energies[2*(order+1)+k_order] += q1*c2*XCOrderD1(λ2,aux_dist);
                    elseif k_order == 3
                        energies[0*(order+1)+k_order] += q1*c2*XCOrder2(λ2,aux_dist);
                        energies[2*(order+1)+k_order] += q1*c2*XCOrderD2(λ2,aux_dist);
                    elseif k_order == 4
                        energies[0*(order+1)+k_order] += q1*c2*XCOrder3(λ2,aux_dist);
                        energies[2*(order+1)+k_order] += q1*c2*XCOrderD3(λ2,aux_dist);
                    elseif k_order == 5
                        energies[0*(order+1)+k_order] += q1*c2*XCOrder4(λ2,aux_dist);
                        energies[2*(order+1)+k_order] += q1*c2*XCOrderD4(λ2,aux_dist);
                    elseif k_order == 6
                        energies[0*(order+1)+k_order] += q1*c2*XCOrder5(λ2,aux_dist);
                        energies[2*(order+1)+k_order] += q1*c2*XCOrderD5(λ2,aux_dist);
                    elseif k_order == 7
                        energies[0*(order+1)+k_order] += q1*c2*XCOrder6(λ2,aux_dist);
                        energies[2*(order+1)+k_order] += q1*c2*XCOrderD6(λ2,aux_dist);
                    elseif k_order == 8
                        energies[0*(order+1)+k_order] += q1*c2*XCOrder7(λ2,aux_dist);
                        energies[2*(order+1)+k_order] += q1*c2*XCOrderD7(λ2,aux_dist);
                    elseif k_order == 9
                        energies[0*(order+1)+k_order] += q1*c2*XCOrder8(λ2,aux_dist);
                        energies[2*(order+1)+k_order] += q1*c2*XCOrderD8(λ2,aux_dist);
                    elseif k_order == 10
                        energies[0*(order+1)+k_order] += q1*c2*XCOrder9(λ2,aux_dist);
                        energies[2*(order+1)+k_order] += q1*c2*XCOrderD9(λ2,aux_dist);
                    elseif k_order == 11
                        energies[0*(order+1)+k_order] += q1*c2*XCOrder10(λ2,aux_dist);
                        energies[2*(order+1)+k_order] += q1*c2*XCOrderD10(λ2,aux_dist);
                    end
                end
            end
        end

        # cloud-cloud
        for i in 1:num_clouds1
            for j in i:num_clouds1
                c1 = molecule1.cloud_data[i,4];
                λ1 = molecule1.cloud_data[i,5];
                c2 = molecule1.cloud_data[j,4];
                λ2 = molecule1.cloud_data[j,5];

                if (c1 == 0) || (c2 == 0)
                    continue;
                end

                λ = (λ1*λ2)/(λ1+λ2);

                # polarization
                c1 *= molecule1.cloud_data[i,6];
                c2 *= molecule1.cloud_data[j,6];

                aux_dist = molecule1.cloud_data[i,1:3];
                aux_dist -= molecule1.cloud_data[j,1:3];
                aux_dist = norm(aux_dist);

                C1 = copy(c1);
                c1 *= exp(-λ*(aux_dist^2.0));

                aux_dist = zeros(aux_type,22) .+ aux_dist;
                aux_dist = aux_dist .^ collect(1:22);

                λ = zeros(aux_type,22) .+ λ;
                λ = λ .^ collect(1:22);
                for k_order in 1:(order+1)
                    if k_order == 1
                        energies[1*(order+1)+k_order] += C1*c2*XCOrder0(λ,aux_dist);
                        energies[3*(order+1)+k_order] += C1*c2*XCOrderD01(λ,aux_dist);
                        energies[3*(order+1)+k_order] += c1*c2*XCOrderD00(λ,aux_dist);
                    elseif k_order == 2
                        energies[1*(order+1)+k_order] += c1*c2*XCOrder1(λ);
                        energies[3*(order+1)+k_order] += c1*c2*XCOrderD1(λ,aux_dist);
                    elseif k_order == 3
                        energies[1*(order+1)+k_order] += c1*c2*XCOrder2(λ,aux_dist);
                        energies[3*(order+1)+k_order] += c1*c2*XCOrderD2(λ,aux_dist);
                    elseif k_order == 4
                        energies[1*(order+1)+k_order] += c1*c2*XCOrder3(λ,aux_dist);
                        energies[3*(order+1)+k_order] += c1*c2*XCOrderD3(λ,aux_dist);
                    elseif k_order == 5
                        energies[1*(order+1)+k_order] += c1*c2*XCOrder4(λ,aux_dist);
                        energies[3*(order+1)+k_order] += c1*c2*XCOrderD4(λ,aux_dist);
                    elseif k_order == 6
                        energies[1*(order+1)+k_order] += c1*c2*XCOrder5(λ,aux_dist);
                        energies[3*(order+1)+k_order] += c1*c2*XCOrderD5(λ,aux_dist);
                    elseif k_order == 7
                        energies[1*(order+1)+k_order] += c1*c2*XCOrder6(λ,aux_dist);
                        energies[3*(order+1)+k_order] += c1*c2*XCOrderD6(λ,aux_dist);
                    elseif k_order == 8
                        energies[1*(order+1)+k_order] += c1*c2*XCOrder7(λ,aux_dist);
                        energies[3*(order+1)+k_order] += c1*c2*XCOrderD7(λ,aux_dist);
                    elseif k_order == 9
                        energies[1*(order+1)+k_order] += c1*c2*XCOrder8(λ,aux_dist);
                        energies[3*(order+1)+k_order] += c1*c2*XCOrderD8(λ,aux_dist);
                    elseif k_order == 10
                        energies[1*(order+1)+k_order] += c1*c2*XCOrder9(λ,aux_dist);
                        energies[3*(order+1)+k_order] += c1*c2*XCOrderD9(λ,aux_dist);
                    elseif k_order == 11
                        energies[1*(order+1)+k_order] += c1*c2*XCOrder10(λ,aux_dist);
                        energies[3*(order+1)+k_order] += c1*c2*XCOrderD10(λ,aux_dist);
                    end
                end
            end
        end

        # Intermolecular interactions
        for jj in (ii+1):num_molecules
            molecule2 = molecules[jj];
            num_atoms2 = size(molecule2.atoms_data)[1];
            num_clouds2 = size(molecule2.cloud_data)[1];

            # nuclei-cloud 
            for i in 1:num_atoms1
                for j in 1:num_clouds2
                    q1 = molecule1.atoms_data[i,4];
                    c2 = molecule2.cloud_data[j,4];
                    λ2 = molecule2.cloud_data[j,5];

                    if c2 == 0
                        continue;
                    end

                    # polarization
                    c2 *= molecule2.cloud_data[j,6];

                    aux_dist = molecule1.atoms_data[i,1:3];
                    aux_dist -= molecule2.cloud_data[j,1:3];
                    aux_dist = norm(aux_dist);

                    C2 = copy(c2);
                    c2 *= exp(-λ2*(aux_dist^2.0));

                    aux_dist = zeros(aux_type,22) .+ aux_dist;
                    aux_dist = aux_dist .^ collect(1:22);

                    λ2 = zeros(aux_type,22) .+ λ2;
                    λ2 = λ2 .^ collect(1:22);
                    for k_order in 1:(order+1)
                        if k_order == 1
                            energies[0*(order+1)+k_order] += q1*C2*XCOrder0(λ2,aux_dist);
                            energies[2*(order+1)+k_order] += q1*C2*XCOrderD01(λ2,aux_dist);
                            energies[2*(order+1)+k_order] += q1*c2*XCOrderD00(λ2,aux_dist);
                        elseif k_order == 2
                            energies[0*(order+1)+k_order] += q1*c2*XCOrder1(λ2);
                            energies[2*(order+1)+k_order] += q1*c2*XCOrderD1(λ2,aux_dist);
                        elseif k_order == 3
                            energies[0*(order+1)+k_order] += q1*c2*XCOrder2(λ2,aux_dist);
                            energies[2*(order+1)+k_order] += q1*c2*XCOrderD2(λ2,aux_dist);
                        elseif k_order == 4
                            energies[0*(order+1)+k_order] += q1*c2*XCOrder3(λ2,aux_dist);
                            energies[2*(order+1)+k_order] += q1*c2*XCOrderD3(λ2,aux_dist);
                        elseif k_order == 5
                            energies[0*(order+1)+k_order] += q1*c2*XCOrder4(λ2,aux_dist);
                            energies[2*(order+1)+k_order] += q1*c2*XCOrderD4(λ2,aux_dist);
                        elseif k_order == 6
                            energies[0*(order+1)+k_order] += q1*c2*XCOrder5(λ2,aux_dist);
                            energies[2*(order+1)+k_order] += q1*c2*XCOrderD5(λ2,aux_dist);
                        elseif k_order == 7
                            energies[0*(order+1)+k_order] += q1*c2*XCOrder6(λ2,aux_dist);
                            energies[2*(order+1)+k_order] += q1*c2*XCOrderD6(λ2,aux_dist);
                        elseif k_order == 8
                            energies[0*(order+1)+k_order] += q1*c2*XCOrder7(λ2,aux_dist);
                            energies[2*(order+1)+k_order] += q1*c2*XCOrderD7(λ2,aux_dist);
                        elseif k_order == 9
                            energies[0*(order+1)+k_order] += q1*c2*XCOrder8(λ2,aux_dist);
                            energies[2*(order+1)+k_order] += q1*c2*XCOrderD8(λ2,aux_dist);
                        elseif k_order == 10
                            energies[0*(order+1)+k_order] += q1*c2*XCOrder9(λ2,aux_dist);
                            energies[2*(order+1)+k_order] += q1*c2*XCOrderD9(λ2,aux_dist);
                        elseif k_order == 11
                            energies[0*(order+1)+k_order] += q1*c2*XCOrder10(λ2,aux_dist);
                            energies[2*(order+1)+k_order] += q1*c2*XCOrderD10(λ2,aux_dist);
                        end
                    end
                end
            end

            for i in 1:num_atoms2
                for j in 1:num_clouds1
                    c1 = molecule1.cloud_data[j,4];
                    λ1 = molecule1.cloud_data[j,5];
                    q2 = molecule2.atoms_data[i,4];

                    if c1 == 0
                        continue;
                    end

                    # polarization
                    c1 *= molecule1.cloud_data[j,6];

                    aux_dist = molecule2.atoms_data[i,1:3];
                    aux_dist -= molecule1.cloud_data[j,1:3];
                    aux_dist = norm(aux_dist);

                    C1 = copy(c1);
                    c1 *= exp(-λ1*(aux_dist^2.0));

                    aux_dist = zeros(aux_type,22) .+ aux_dist;
                    aux_dist = aux_dist .^ collect(1:22);

                    λ1 = zeros(aux_type,22) .+ λ1;
                    λ1 = λ1 .^ collect(1:22);
                    for k_order in 1:(order+1)
                        if k_order == 1
                            energies[0*(order+1)+k_order] += q2*C1*XCOrder0(λ1,aux_dist);
                            energies[2*(order+1)+k_order] += q2*C1*XCOrderD01(λ1,aux_dist);
                            energies[2*(order+1)+k_order] += q2*c1*XCOrderD00(λ1,aux_dist);
                        elseif k_order == 2
                            energies[0*(order+1)+k_order] += q2*c1*XCOrder1(λ1);
                            energies[2*(order+1)+k_order] += q2*c1*XCOrderD1(λ1,aux_dist);
                        elseif k_order == 3
                            energies[0*(order+1)+k_order] += q2*c1*XCOrder2(λ1,aux_dist);
                            energies[2*(order+1)+k_order] += q2*c1*XCOrderD2(λ1,aux_dist);
                        elseif k_order == 4
                            energies[0*(order+1)+k_order] += q2*c1*XCOrder3(λ1,aux_dist);
                            energies[2*(order+1)+k_order] += q2*c1*XCOrderD3(λ1,aux_dist);
                        elseif k_order == 5
                            energies[0*(order+1)+k_order] += q2*c1*XCOrder4(λ1,aux_dist);
                            energies[2*(order+1)+k_order] += q2*c1*XCOrderD4(λ1,aux_dist);
                        elseif k_order == 6
                            energies[0*(order+1)+k_order] += q2*c1*XCOrder5(λ1,aux_dist);
                            energies[2*(order+1)+k_order] += q2*c1*XCOrderD5(λ1,aux_dist);
                        elseif k_order == 7
                            energies[0*(order+1)+k_order] += q2*c1*XCOrder6(λ1,aux_dist);
                            energies[2*(order+1)+k_order] += q2*c1*XCOrderD6(λ1,aux_dist);
                        elseif k_order == 8
                            energies[0*(order+1)+k_order] += q2*c1*XCOrder7(λ1,aux_dist);
                            energies[2*(order+1)+k_order] += q2*c1*XCOrderD7(λ1,aux_dist);
                        elseif k_order == 9
                            energies[0*(order+1)+k_order] += q2*c1*XCOrder8(λ1,aux_dist);
                            energies[2*(order+1)+k_order] += q2*c1*XCOrderD8(λ1,aux_dist);
                        elseif k_order == 10
                            energies[0*(order+1)+k_order] += q2*c1*XCOrder9(λ1,aux_dist);
                            energies[2*(order+1)+k_order] += q2*c1*XCOrderD9(λ1,aux_dist);
                        elseif k_order == 11
                            energies[0*(order+1)+k_order] += q2*c1*XCOrder10(λ1,aux_dist);
                            energies[2*(order+1)+k_order] += q2*c1*XCOrderD10(λ1,aux_dist);
                        end
                    end
                end
            end

            # cloud-cloud
            for i in 1:num_clouds1
                for j in 1:num_clouds2
                    c1 = molecule1.cloud_data[i,4];
                    λ1 = molecule1.cloud_data[i,5];
                    c2 = molecule2.cloud_data[j,4];
                    λ2 = molecule2.cloud_data[j,5];

                    if (c1 == 0) || (c2 == 0)
                        continue;
                    end

                    λ = (λ1*λ2)/(λ1+λ2);

                    # polarization
                    c1 *= molecule1.cloud_data[i,6];
                    c2 *= molecule2.cloud_data[j,6];

                    aux_dist = molecule1.cloud_data[i,1:3];
                    aux_dist -= molecule2.cloud_data[j,1:3];
                    aux_dist = norm(aux_dist);

                    C1 = copy(c1);
                    c1 *= exp(-λ*(aux_dist^2.0));

                    aux_dist = zeros(aux_type,22) .+ aux_dist;
                    aux_dist = aux_dist .^ collect(1:22);

                    λ = zeros(aux_type,22) .+ λ;
                    λ = λ .^ collect(1:22);
                    for k_order in 1:(order+1)
                        if k_order == 1
                            energies[1*(order+1)+k_order] += C1*c2*XCOrder0(λ,aux_dist);
                            energies[3*(order+1)+k_order] += C1*c2*XCOrderD01(λ,aux_dist);
                            energies[3*(order+1)+k_order] += c1*c2*XCOrderD00(λ,aux_dist);
                        elseif k_order == 2
                            energies[1*(order+1)+k_order] += c1*c2*XCOrder1(λ);
                            energies[3*(order+1)+k_order] += c1*c2*XCOrderD1(λ,aux_dist);
                        elseif k_order == 3
                            energies[1*(order+1)+k_order] += c1*c2*XCOrder2(λ,aux_dist);
                            energies[3*(order+1)+k_order] += c1*c2*XCOrderD2(λ,aux_dist);
                        elseif k_order == 4
                            energies[1*(order+1)+k_order] += c1*c2*XCOrder3(λ,aux_dist);
                            energies[3*(order+1)+k_order] += c1*c2*XCOrderD3(λ,aux_dist);
                        elseif k_order == 5
                            energies[1*(order+1)+k_order] += c1*c2*XCOrder4(λ,aux_dist);
                            energies[3*(order+1)+k_order] += c1*c2*XCOrderD4(λ,aux_dist);
                        elseif k_order == 6
                            energies[1*(order+1)+k_order] += c1*c2*XCOrder5(λ,aux_dist);
                            energies[3*(order+1)+k_order] += c1*c2*XCOrderD5(λ,aux_dist);
                        elseif k_order == 7
                            energies[1*(order+1)+k_order] += c1*c2*XCOrder6(λ,aux_dist);
                            energies[3*(order+1)+k_order] += c1*c2*XCOrderD6(λ,aux_dist);
                        elseif k_order == 8
                            energies[1*(order+1)+k_order] += c1*c2*XCOrder7(λ,aux_dist);
                            energies[3*(order+1)+k_order] += c1*c2*XCOrderD7(λ,aux_dist);
                        elseif k_order == 9
                            energies[1*(order+1)+k_order] += c1*c2*XCOrder8(λ,aux_dist);
                            energies[3*(order+1)+k_order] += c1*c2*XCOrderD8(λ,aux_dist);
                        elseif k_order == 10
                            energies[1*(order+1)+k_order] += c1*c2*XCOrder9(λ,aux_dist);
                            energies[3*(order+1)+k_order] += c1*c2*XCOrderD9(λ,aux_dist);
                        elseif k_order == 11
                            energies[1*(order+1)+k_order] += c1*c2*XCOrder10(λ,aux_dist);
                            energies[3*(order+1)+k_order] += c1*c2*XCOrderD10(λ,aux_dist);
                        end
                    end
                end
            end
        end
    end

    return energies;
end
