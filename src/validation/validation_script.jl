using BiochemicalAlgorithms, DataFrames, BenchmarkTools

export export_all_gaff_paper_files_to_mol2, export_all_pdb_test_files_to_mol2, compare_mol_antechamber_to_balljl, load_multiCompound_pubchem_json_and_export_to_mol2, 
    export_all_from_directory_with_gaff, export_all_from_directory, atomtype_comparison, count_assigned_properties, count_assigned_property_fields, 
    load_all_gaff_paper_files, run_atomtyping_on_gaff_paper_files, atomtype_all_gaff_mol2_files, load_all_from_directory, benchmark_all_def_files, print_table_PubChem_ids

function print_table_PubChem_ids(folder::String)
    for filename in sort(readdir(folder))
        println(filename[1:findfirst('_', filename)-1], " & ", filename[findlast('_', filename)+1:findlast('.', filename)-1], " \\\\")
    end
end


function benchmark_all_def_files()
    data_folder = "data/antechamber/ATOMTYPE_"
    def_vec = ["DU", "SYBYL", "GFF", "AMBER"]
    FDA_mol2_path = "../huettel-msc/FDA_PubChem_category_mol2_files/"
    benchmarks_path = string(FDA_mol2_path, "benchmarks/")
    export_all_from_directory("../huettel-msc/extended_test_cases/", "../huettel-msc/extended_test_cases/", string(data_folder, "DU.DEF"))
    for def in def_vec
        println(def)
        deffile = string(data_folder, def, ".DEF")
        export_to_folder = string(benchmarks_path, def, "/")
        bench = @benchmarkable export_all_from_directory($FDA_mol2_path, $export_to_folder, $deffile)
        display(run(bench, samples = 7, seconds = 120))
    end
end


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
            # export_mol2(mol, export_folder)
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
    # for atmNum in result_list 
    #     println(basename(mol1.name), "   atmNum: ", atmNum, "  ", " ball atomtype: ", mol1.atoms.atomtype[atmNum], "  ; ", " julia atomtype: ", mol2.atoms.atomtype[atmNum])
    # end
    return result_list
end


function count_assigned_properties(def_df::DataFrame)
    counter = 0
    for col in eachcol(def_df)[2:end-2]
        counter += count(x -> x != nothing, col)
    end
    for col in eachcol(def_df)[end-1:end]
        for colrow in filter(x -> x != nothing, col)
            counter += count(x -> x in ['[','(','.', ','], colrow)
        end
    end
    return counter
end


function count_assigned_property_fields(def_df::DataFrame)
    counter = 0
    for col in eachcol(def_df)[2:end]
        counter += count(x -> x != nothing, col)
    end
    return counter
end