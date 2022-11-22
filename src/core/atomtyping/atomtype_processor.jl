using BiochemicalAlgorithms

function get_molecule_atomtypes!(mol::AbstractMolecule, mapfile::AbstractString)
    PreprocessingMolecule!(mol)
    def_file = load_atomtyping_DEF(mapfile)
    
    # loop over atoms, prefilter dataframe by element and neighbor count
    for i = (1:nrow(mol.atoms))
        df_curr = copy(def_file)
        EWG_groups = -1
        if Int(mol.atoms.element[i]) == 1
            EWG_groups = count_EWG(i, mol)
        end
        if !isempty(df_curr[(df_curr[!,:atomic_number] .== Int(mol.atoms.element[i])), :])
            df_curr = df_curr[(df_curr[!,:atomic_number] .== Int(mol.atoms.element[i])), :]
        end
        if !isempty(df_curr[(df_curr[!,:num_neighbors] .== lastindex(neighbors(mol.properties["mol_graph"], i))), :])
            df_curr = df_curr[(df_curr[!,:num_neighbors] .== lastindex(neighbors(mol.properties["mol_graph"], i))), :]
        end
        if !isempty(df_curr[(df_curr[!,:num_H_bonds] .== neighbors(mol.properties[i]["mol_graph"], i)), :])

        end
        df_curr = df_curr[(df_curr[!,:atomic_number] .== Int(mol.atoms.element[i])) .& (df_curr[!,:num_neighbors] .== Int(mol.atoms.num_neighbors[i])),:]

    end
end


function count_EWG(num::Int, mol::AbstractMolecule)
    # Electron withdrawal Atoms according to antechamber document are only: N, O, F, Cl, and Br
    # To Do: Test differences for Atom Typing, see below typical know EWG
    strong_pullers = ["Cl1", "F1", "Br1", "I1", "O1", "S1"]
    possible_pullers = ["C2", "C3", "C4", "S3", "N3", "P3", "P4", "O2", "S2"]
    elec_pullers_num = 0
    if true in in(mol.atoms.properties[num]["ElementWithNeighborCount"]).(possible_pullers)
        for i = (1:lastindex(neighbors(mol.properties["mol_graph"], num)))
            if true in in(strong_pullers).(ATD_df.Element_wNeighborCount[ATD_df.Secondary_Neighbors[num][i]])
                # all neighbors of neighbor that are in strong_pullers are an EWG
                elec_pullers_num += countmap(in(strong_pullers).(mol.atoms.elem[ATD_df.Secondary_Neighbors[num][i]]))[true]
            end
        end
    end
    if elec_pullers_num == 0
        elec_pullers_num = -1
    end
    # # println("EWG: ", elec_pullers_num)
    return elec_pullers_num
end