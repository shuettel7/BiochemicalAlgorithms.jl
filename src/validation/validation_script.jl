using BiochemicalAlgorithms, DataFrames

export export_all_gaff_paper_files_to_mol2, export_all_pdb_test_files_to_mol2, compare_mol_antechamber_to_balljl


function export_all_gaff_paper_files_to_mol2()
    mol_df = load_all_gaff_paper_files()
    for num = (1:nrow(mol_df))
        gaff_atomtyping_wrapper!(mol_df.abstract_mol[num])
        export_mol2(mol_df.abstract_mol[num], "../huettel-msc/export_folder/validation_script_testing_folder/")
    end
end


function export_all_pdb_test_files_to_mol2(file_location_a::AbstractString, export_to_dir_location::AbstractString)
    mol_df = load_all_pdb_test_files(file_location_a)
    for num = (1:nrow(mol_df))
        export_mol2(mol_df.abstract_mol[num], export_to_dir_location)
    end
end

function compare_mol_antechamber_to_balljl(directory1::AbstractString, directory2::AbstractString)
    comparison_df = DataFrame([Vector{String}(), Vector{Vector{Any}}()], ["molname", "values"])
    for i = (1:lastindex(readdir(directory1)))
        filename = basename(readdir(directory1)[i])
        mol1 = load_mol2(string(directory1, filename))
        mol2 = load_mol2(string(directory2, filename))
        atomtype_comparison_vector = atomtype_comparison(mol1, mol2)
        push!(comparison_df, (basename(mol1.name), atomtype_comparison_vector))
    end
    return comparison_df
end

function atomtype_comparison(mol1::AbstractMolecule, mol2::AbstractMolecule)
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