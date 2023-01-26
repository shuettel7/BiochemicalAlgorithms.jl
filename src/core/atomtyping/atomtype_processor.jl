using DataFramesMeta, DataFrames

export get_molecule_atomtypes!, gaff_atomtyping_wrapper!

function gaff_atomtyping_wrapper!(mol::Molecule, mapfile = "data/antechamber/ATOMTYPE_GFF.DEF")
    PreprocessingMolecule!(mol)
    get_molecule_atomtypes!(mol, mapfile)
    gaff_postprocessing_all_conjugated_systems!(mol)
end

function get_molecule_atomtypes!(mol::Molecule, mapfile::String)
    def_file_df = load_atomtyping_DEF(mapfile)
    atmprops_df = mol.properties["atmprops_df"]    
    
    # loop over atoms, prefilter dataframe by element and neighbor count
    Threads.@threads for i = (1:nrow(mol.atoms))
        num_H_neighbors = countmap(in(["H1"]).(atmprops_df.ElementWithNeighborCount[atmprops_df.Neighbors[i]]), alg=:dict)[true]
        num_EWG_groups = mol.atoms[i, :properties]["num_EWG_groups"]

        for row in eachrow(def_file_df)
            if (isnothing(row.atomic_number) || row.atomic_number == Int(mol.atoms.element[i])) &&
                (isnothing(row.num_neighbors) || row.num_neighbors == lastindex(atmprops_df.Neighbors[i])) &&
                (isnothing(row.num_H_bonds) || row.num_H_bonds == num_H_neighbors) &&
                (isnothing(row.electron_withdrawal_groups) || row.electron_withdrawal_groups == num_EWG_groups) &&
                (isnothing(row.atomic_property) || APS_processor(row.atomic_property, atmprops_df[i,:])) &&
                (isnothing(row.CES) || CES_processor(CES_parser(row.CES, mol, i), mol, i))
                mol.atoms.atomtype[i] = row.type_name
                break
            end
        end
    end
end