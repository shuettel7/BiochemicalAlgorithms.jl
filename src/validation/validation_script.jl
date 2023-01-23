using BiochemicalAlgorithms, DataFrames

export export_all_gaff_paper_files_to_mol2, export_all_pdb_test_files_to_mol2, compare_mol_antechamber_to_balljl, load_multiCompound_pubchem_json_and_export_to_mol2, 
    export_all_from_directory_with_gaff, export_all_from_directory


function export_all_from_directory_with_gaff(directory::String, toDirectory::String)
    mol_df = load_all_from_directory(directory)
    for num = (1:nrow(mol_df))
        gaff_atomtyping_wrapper!(mol_df.abstract_mol[num])
        export_mol2(mol_df.abstract_mol[num], toDirectory)
    end
end

function export_all_from_directory(directory::String, toDirectory::String, def_file::String)
    mol_df = load_all_from_directory(directory)
    if contains(def_file, "ATOMTYPE_GFF")
        Threads.@threads for num = (1:nrow(mol_df))
            PreprocessingMolecule!(mol_df.abstract_mol[num])
            get_molecule_atomtypes!(mol_df.abstract_mol[num], def_file)
            gaff_postprocessing_all_conjugated_systems!(mol_df.abstract_mol[num])
            export_mol2(mol_df.abstract_mol[num], toDirectory)
        end
    else
        Threads.@threads for num =  (1:nrow(mol_df))
            PreprocessingMolecule!(mol_df.abstract_mol[num])
            get_molecule_atomtypes!(mol_df.abstract_mol[num], def_file)
            export_mol2(mol_df.abstract_mol[num], toDirectory)
        end
    end
end

function load_multiCompound_pubchem_json_and_export_to_mol2(fname::String, export_folder::String)
    for (i,line) in enumerate(readlines(fname))
        if i != 1 && i != lastindex(readlines(fname))
            mol = load_pubchem_json(line)
            export_mol2(mol, export_folder)
        end
    end
end


function export_all_gaff_paper_files_to_mol2(toDirectory::String)
    mol_df = load_all_gaff_paper_files()
    for num = (1:nrow(mol_df))
        gaff_atomtyping_wrapper!(mol_df.abstract_mol[num])
        export_mol2(mol_df.abstract_mol[num], toDirectory)
    end
end


function export_all_pdb_test_files_to_mol2(file_location_a::String, export_to_dir_location::String)
    mol_df = load_all_pdb_test_files(file_location_a)
    for num = (1:nrow(mol_df))
        export_mol2(mol_df.abstract_mol[num], export_to_dir_location)
    end
end

function compare_mol_antechamber_to_balljl(directory1::String, directory2::String)
    comparison_df = DataFrame([Vector{String}(), Vector{Vector{Any}}()], ["molname", "values"])
    for i = (1:lastindex(readdir(directory1)))
        filename = basename(readdir(directory1)[i])
        mol1 = load_mol2(string(directory1, filename))
        mol2 = load_mol2(string(directory2, filename))
        atomtype_comparison_vector = (isnothing(mol1) || isnothing(mol2)) ? ["missing"] : atomtype_comparison(mol1, mol2)
        push!(comparison_df, (isnothing(mol1) ? basename(mol2.name) : basename(mol1.name), atomtype_comparison_vector))
    end
    return comparison_df
end

function atomtype_comparison(mol1::Molecule, mol2::Molecule)
    result_list = Vector{Any}()
    if nrow(mol1.atoms) != nrow(mol2.atoms)
        return "Molecules are of different size"
    else
        for i = (1:nrow(mol1.atoms))
            if string(mol1.atoms.atomtype[i]) != string(mol2.atoms.atomtype[i])
                append!(result_list, i)
            end
        end
    end
    if isempty(result_list)
        return [true]
    end
    return result_list
end