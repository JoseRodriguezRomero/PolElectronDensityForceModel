using Printf

include("GetPotentialFromAngles.jl")

function MinimumDistanceMolecs(mol1::Molecule, mol2::Molecule)
    # Returs the distance between the closest two atoms between 
    # mol1 and mol2.
    min_dist = mol2.atoms_data[1,1:3] - mol1.atoms_data[1,1:3];
    min_dist = sqrt(sum(min_dist.^2.0));

    num_atoms1 = size(mol1.atoms_data)[1];
    num_atoms2 = size(mol2.atoms_data)[1];

    for i in 1:num_atoms1
        for j in 1:num_atoms2
            r1 = mol1.atoms_data[i,1:3];
            r2 = mol2.atoms_data[j,1:3];

            aux_dist = sqrt(sum((r1-r2).^2.0));
            min_dist = minimum([min_dist,aux_dist]);
        end
    end
    return min_dist;
end

function MakeData(base_dir::String, num_samples::Int)
    # Makes num_samples number of fitting data by rotating and 
    # translating a molecule of water (one stays put).
    density_file = "H2O_ecp_fitted_data.txt";

    fileID1 = open(base_dir*"RotationsAndDisplacements.txt","w");
    mol_a = ReadMolecule(density_file);
    CenterAtAtomIndex!(mol_a,1); # Position the Oxygen atom at (0,0,0)

    @printf fileID1 "%18s %18s %18s " "theta_x" "theta_y" "theta_z";
    @printf fileID1 "%18s %18s %18s \n" "dx" "dy" "dz";

    for i in 1:num_samples
        # OH Coord
        ang_samples = 25;
        # for j in 1:ang_samples
        #     mol_b = copy(mol_a);
        #     angs_and_disp = zeros(Float64,1,6);
        #     θ = - (π/2.0) * ((j-1)/(ang_samples-1));

        #     angs_and_disp[3] = θ;
        #     angs_and_disp[4:6] += mol_a.atoms_data[2,1:3];
        #     angs_and_disp[4:6] -= mol_a.atoms_data[1,1:3];
        #     OH_dist = sqrt(sum(angs_and_disp[4:6].^2.0));
        #     angs_and_disp[4:6] ./= OH_dist;
        #     angs_and_disp[4:6] .*= OH_dist + 1.0 + (i-1)/(num_samples-1);

        #     MoveAndRotateMolec!(mol_b,angs_and_disp);

        #     for k in 1:6
        #         @printf fileID1 "%18.12f " angs_and_disp[k];
        #     end
                
        #     @printf fileID1 "\n";
        # end

        # OO Coord
        ang_samples = 25;
        for j in 1:ang_samples
            angs_and_disp = zeros(Float64,1,6);
        
            mol_b = copy(mol_a);
            θ = - (π/2.0) * ((j-1)/(ang_samples-1));
            angs_and_disp[3] = π - θ;
            angs_and_disp[5] = 2.0 + 10.0*(i-1)/(num_samples-1);
            MoveAndRotateMolec!(mol_b,angs_and_disp);
        
            for k in 1:6
                @printf fileID1 "%18.12f " angs_and_disp[k];
            end
        
            @printf fileID1 "\n";
        end
    end

    close(fileID1);
end

MakeData("Training Data/",200)
