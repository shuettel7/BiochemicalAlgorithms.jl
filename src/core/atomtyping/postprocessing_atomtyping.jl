
export gaff_postprocessing_all_conjugated_systems!

function gaff_postprocessing_all_conjugated_systems!(mol::Molecule)
    mol_graph = mol.properties["mol_graph"]
    group_atomtypes_dict = Dict{Int, Vector{String}}(
                                1 => ["cc", "ce", "cp", "nc", "ne", "pc", "pe"],
                                2 => ["cd", "cf", "cq", "nd", "nf", "pd", "pf"])

    group1_conjugated_atoms = mol.atoms[(in(group_atomtypes_dict[1]).(mol.atoms[!,:atomtype])),:]

    seen_atoms_vec = Vector{Int}()
    
    for conjAtm in eachrow(group1_conjugated_atoms)
        if conjAtm.number in seen_atoms_vec
            continue
        end
        neighbors_in_conjugated_atoms = mol.atoms[(in(group_atomtypes_dict[1]).(mol.atoms[!,:atomtype]) .&& 
                                            in(neighbors(mol_graph, conjAtm.number)).(mol.atoms[!,:number])),:]
        
        depth = 1
        next_neigh_vec = neighbors_in_conjugated_atoms.number
        start_atom = conjAtm.number
        for neigh in next_neigh_vec
            if neigh in seen_atoms_vec
                continue
            end
            all_paths_from_curr_conjAtm = conjAtm_path_builder!(group1_conjugated_atoms, [start_atom, neigh], mol, depth)
            append!(seen_atoms_vec, reduce(vcat, all_paths_from_curr_conjAtm))
        end
    end
end


function conjAtm_path_builder!(group1_conjugated_atoms::DataFrame, path_vec::Vector{Int}, mol::Molecule, depth::Int64)
    mol_graph = mol.properties["mol_graph"]
    group_atomtypes_dict = Dict{Int, Vector{String}}(
                                1 => ["cc", "ce", "cp", "nc", "ne", "pc", "pe"],
                                2 => ["cd", "cf", "cq", "nd", "nf", "pd", "pf"])
    source_atom = path_vec[1]
    previous_atom = path_vec[lastindex(path_vec)-1]
    curr_atom = path_vec[lastindex(path_vec)]

    combination_df = DataFrame(:a1 => [previous_atom, curr_atom], :a2 => [curr_atom, previous_atom])
    bond_in_conjugated_atoms = innerjoin(combination_df, mol.bonds, on = [:a1, :a2])

    if haskey(bond_in_conjugated_atoms.properties[1], "TRIPOS_tag") && bond_in_conjugated_atoms.properties[1]["TRIPOS_tag"] == "ar"
        if mol.atoms.atomtype[previous_atom] in group_atomtypes_dict[1] && in(mol.atoms.atomtype[curr_atom], group_atomtypes_dict[1])
            mol.atoms.atomtype[curr_atom] = group_atomtypes_dict[2][findfirst(x -> x == mol.atoms.atomtype[curr_atom], group_atomtypes_dict[1])]
        end
    elseif bond_in_conjugated_atoms.order[1] == BondOrder.T(1)
        if mol.atoms.atomtype[previous_atom] in group_atomtypes_dict[2] && in(mol.atoms.atomtype[curr_atom], group_atomtypes_dict[1])
            mol.atoms.atomtype[curr_atom] = group_atomtypes_dict[2][findfirst(x -> x == mol.atoms.atomtype[curr_atom], group_atomtypes_dict[1])]
        end
    elseif bond_in_conjugated_atoms.order[1] == BondOrder.T(2)
        if mol.atoms.atomtype[previous_atom] in group_atomtypes_dict[1] && in(mol.atoms.atomtype[curr_atom], group_atomtypes_dict[1])
            mol.atoms.atomtype[curr_atom] = group_atomtypes_dict[2][findfirst(x -> x == mol.atoms.atomtype[curr_atom], group_atomtypes_dict[1])]
        end
    end
    next_neighbors_vec = filter(x -> !in(x, path_vec) && x in neighbors(mol_graph, curr_atom) && 
                                    x in group1_conjugated_atoms.number, neighborhood(mol_graph, source_atom, depth+1))
                    
    all_paths = Vector{Vector{Int}}()
    if isempty(next_neighbors_vec)
        push!(all_paths, path_vec)
        return all_paths
    else
        for nextNeigh in next_neighbors_vec
            append!(all_paths, conjAtm_path_builder!(group1_conjugated_atoms, vcat(path_vec, nextNeigh), mol, depth+1))
        end    
    end
    return all_paths     
end