using DataFramesMeta, DataFrames

export get_molecule_atomtypes!, gaff_atomtyping_wrapper!

function gaff_atomtyping_wrapper!(mol::AbstractMolecule, mapfile = "data/antechamber/ATOMTYPE_GFF.DEF")
    PreprocessingMolecule!(mol)
    get_molecule_atomtypes!(mol, mapfile)
    gaff_postprocessing_all_conjugated_systems!(mol)
end

function get_molecule_atomtypes!(mol::AbstractMolecule, mapfile::AbstractString)
    def_file_df = load_atomtyping_DEF(mapfile)
    atmprops_df = mol.properties["atmprops_df"]
    
    # loop over atoms, prefilter dataframe by element and neighbor count
    for i = (1:nrow(mol.atoms))
        def_curr_df = copy(def_file_df)
        num_H_neighbors = countmap(in(["H1"]).(atmprops_df.ElementWithNeighborCount[atmprops_df.Neighbors[i]]))[true]
        num_EWG_groups = mol.properties["atmprops_df"].num_EWG_groups[i]
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
                    if colnum == 6 && CES_processor(CES_parser(coldata, mol, i), atmprops_df[i,:], mol, i)
                        match_list[colnum] = 1
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