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

        # start with prefiltering for more easily checkable properties
        def_curr_df = @subset def_curr_df @byrow begin
            (isnothing(:residue_names) || (typeof(mol) == PDBMolecule && :residue_names == mol.atoms.residue_name[i]))
            (isnothing(:atomic_number) || :atomic_number == Int(mol.atoms.element[i])) &&
            (isnothing(:num_neighbors) || :num_neighbors == lastindex(atmprops_df.Neighbors[i])) &&
            (isnothing(:num_H_bonds) || :num_H_bonds == num_H_neighbors) &&
            (isnothing(:electron_withdrawal_groups) || :electron_withdrawal_groups == num_EWG_groups) 
        end

        # only do atom_property and CES checks if there are more atomtypes 
        # than the "unknown atomtype" ("DU" for GAFF) in def_curr_df
        if nrow(def_curr_df) > 1 
            def_curr_df = @subset def_curr_df @byrow begin
                (isnothing(:atomic_property) || APS_processor(:atomic_property, atmprops_df[i,:])) &&
                (isnothing(:CES) || CES_processor(CES_parser(:CES, mol, i), atmprops_df, mol, i))
            end
        end
        mol.atoms.atomtype[i] = def_curr_df.type_name[1]
    end
end