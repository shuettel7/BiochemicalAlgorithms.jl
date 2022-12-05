using BiochemicalAlgorithms, DataFramesMeta, DataFrames

export get_molecule_atomtypes!, count_EWG

function get_molecule_atomtypes!(mol::AbstractMolecule, mapfile::AbstractString)
    PreprocessingMolecule!(mol)
    atmprops_df = create_atom_preprocessing_df!(mol)
    def_file_df = load_atomtyping_DEF(mapfile)
    
    # loop over atoms, prefilter dataframe by element and neighbor count
    for i = (1:nrow(mol.atoms))
        def_curr_df = copy(def_file_df)
        num_H_neighbors = countmap(in(["H1"]).(atmprops_df.ElementWithNeighborCount[atmprops_df.Neighbors[i]]))[true]
        num_EWG_groups = -1
        if Int(mol.atoms.element[i]) == 1
            num_EWG_groups = count_EWG(i, mol, atmprops_df)
        end
        def_curr_df = @subset def_curr_df @byrow begin
            :atomic_number == Int(mol.atoms.element[i])
        end
        df_match = false
        while !df_match && nrow(def_curr_df) > 0
            match_list = [0 for _ in 1:6]
            for (colnum, coldata) in enumerate(eachrow(def_curr_df)[1][3:8])
                if coldata == "*" || coldata == -1
                    match_list[colnum] = 2
                    continue
                else
                    if colnum == 1 && coldata == Int(mol.atoms.element[i])
                        match_list[colnum] = 1
                    end
                    if colnum == 2 && coldata == lastindex(atmprops_df.Neighbors[i])
                        match_list[colnum] = 1
                    end
                    if colnum == 3 && coldata == num_H_neighbors
                        match_list[colnum] = 1
                    end
                    if colnum == 4 && coldata == num_EWG_groups
                        match_list[colnum] = 1
                    end
                    if colnum == 5 && APS_processor(coldata, atmprops_df[i,:])
                        match_list[colnum] = 1
                    end
                    if colnum == 6 && CES_parser(coldata, atmprops_df[i,:], mol)
                        match_list[colnum] = 2
                    end
                end
            end
            if all(in([1, 2]).(match_list))
                df_match = true
                mol.atoms.atomtype[i] = def_curr_df.type_name[1]
            elseif nrow(def_curr_df) == 0
                mol.atoms.atomtype[i] = "DU" # DU is from TRIPOS standard. maybe different expression?
            else
                def_curr_df = def_curr_df[2:nrow(def_curr_df), :]
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