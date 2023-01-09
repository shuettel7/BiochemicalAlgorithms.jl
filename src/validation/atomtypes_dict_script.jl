using BiochemicalAlgorithms

export load_all_gaff_paper_files, run_atomtyping_on_gaff_paper_files, atomtype_all_gaff_mol2_files

function load_all_gaff_paper_files()
    file_location_a = "../huettel-msc/exchange_save_folder/data/gaff_paper_examples/a/"
    file_location_b = "../huettel-msc/exchange_save_folder/data/gaff_paper_examples/b/"
    file_location_c = "../huettel-msc/exchange_save_folder/data/gaff_paper_examples/c/"
    row_count_a = lastindex(readdir(file_location_a))
    row_count_b = lastindex(readdir(file_location_b))
    row_count_total = lastindex(readdir(file_location_a)) + lastindex(readdir(file_location_b)) + lastindex(readdir(file_location_c))
    mol_df = DataFrame([Vector{String}(undef, row_count_total), Vector{AbstractMolecule}(undef, row_count_total)], ["molname", "abstract_mol"])
    for (num, i) in enumerate(readdir(file_location_a))
        mol_df.molname[num] = string("mol_a_", i[1:2])
        mol_df.abstract_mol[num] = load_pubchem_json(string(file_location_a, i))
    end
    for (num, i) in enumerate(readdir(file_location_b))
        num_b = num + row_count_a
        mol_df.molname[num_b] = string("mol_b_", i[1:2])
        mol_df.abstract_mol[num_b] = load_pubchem_json(string(file_location_b, i))
    end
    for (num, i) in enumerate(readdir(file_location_c))
        num_c = num + row_count_a + row_count_b
        mol_df.molname[num_c] = string("mol_c_", i[1:2])
        mol_df.abstract_mol[num_c] = load_pubchem_json(string(file_location_c, i))
    end
    return mol_df
end


function load_all_from_directory(directory::String)
    mol_df = DataFrame([Vector{String}(), Vector{AbstractMolecule}()], ["molname", "abstract_mol"])
    for i in filter(x -> x[end-4:end] == ".mol2", readdir(string(directory)))
        push!(mol_df, (string("mol_a_", i[1:2]), load_mol2(string(directory, i))))
    end
    return mol_df
end