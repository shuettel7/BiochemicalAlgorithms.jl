using BiochemicalAlgorithms

export get_molecule_atomtypes!, count_EWG

function get_molecule_atomtypes!(mol::AbstractMolecule, mapfile::AbstractString)
    PreprocessingMolecule!(mol)
    atmprops_df = create_atom_preprocessing_df(mol)
    def_file_df = load_atomtyping_DEF(mapfile)
    
    # loop over atoms, prefilter dataframe by element and neighbor count
    for i = (1:nrow(mol.atoms))
        def_curr_df = copy(def_file_df)
        default_num_H_neighbors = -1
        num_H_neighbors = countmap(in(["H1"]).(atmprops_df.ElementWithNeighborCount[atmprops_df.Neighbors[i]]))[true]
        num_EWG_groups = default_num_EWG_groups = -1
        if Int(mol.atoms.element[i]) == 1
            EWG_groups = count_EWG(i, mol, atmprops_df)
        end
        ### filter DEF file dataframe by atom properties from preprocessing of molecule
        # filter current dataframe def_curr_df by element number
        if !isempty(def_curr_df[(def_curr_df[!,:atomic_number] .== Int(mol.atoms.element[i])), :])
            df_curr = def_curr_df[(def_curr_df[!,:atomic_number] .== Int(mol.atoms.element[i])), :]
        end

        # filter current dataframe def_curr_df by number of neighbors
        if !isempty(df_curr[(df_curr[!,:num_neighbors] .== lastindex(atmprops_df.Neighbors[i])), :])
            df_curr = df_curr[(df_curr[!,:num_neighbors] .== lastindex(atmprops_df.Neighbors[i])), :]
        end

        # filter current dataframe def_curr_df by number of hydrogens in direct neighborhood
        if !isempty(df_curr[(df_curr[!,:num_H_bonds] .== num_H_neighbors), :])
            df_curr = df_curr[(df_curr[!,:num_H_bonds] .== num_H_neighbors), :]
        else
            if !isempty(df_curr[(df_curr[!,:num_H_bonds] .== default_num_H_neighbors), :])
                df_curr = df_curr[(df_curr[!,:num_H_bonds] .== default_num_H_neighbors), :]
            end
        end

        # filter current dataframe def_curr_df by number of electron-withdrawal groups
        if !isempty(df_curr[(df_curr[!,:electron_withdrawal_groups] .== num_EWG_groups), :])
            df_curr = df_curr[(df_curr[!,:electron_withdrawal_groups] .== num_EWG_groups), :]
        else
            if !isempty(df_curr[(df_curr[!,:electron_withdrawal_groups] .== default_num_EWG_groups), :])
                df_curr = df_curr[(df_curr[!,:electron_withdrawal_groups] .== default_num_EWG_groups), :]
            end
        end
    end
end


function count_EWG(num::Int, mol::AbstractMolecule, atmprops_df::DataFrame)
    # usually only used on hydrogen atoms:
    # Electron withdrawal Atoms according to antechamber document are: N, O, F, Cl, and Br
    # which are bound to the immediate neighbour
    # To Do: Test differences for Atom Typing, see below typical know EWG
    strong_pullers = ["Cl1", "F1", "Br1", "I1", "O1", "S1"]
    possible_pullers = ["C2", "C3", "C4", "S3", "N3", "P3", "P4", "O2", "S2"]
    elec_pullers_num = 0
    if true in in(atmprops_df.ElementWithNeighborCount[atmprops_df.Neighbors[num]]).(possible_pullers)
        for (i, neigh) in enumerate(atmprops_df.Neighbors[num])
            if true in in(strong_pullers).(atmprops_df.ElementWithNeighborCount[atmprops_df.SecondaryNeighbors[num][i]]) &&
                atmprops_df.ElementWithNeighborCount[neigh] in possible_pullers
                # all neighbors of neighbor that are in strong_pullers are an EWG
                elec_pullers_num += countmap(in(strong_pullers).(atmprops_df.ElementWithNeighborCount[atmprops_df.SecondaryNeighbors[num][i]]))[true]
            end
        end
    end
    if elec_pullers_num == 0
        elec_pullers_num = -1
    end
    return elec_pullers_num
end