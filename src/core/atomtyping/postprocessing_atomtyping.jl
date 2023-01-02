
export gaff_postprocessing_all_conjugated_systems!

function gaff_postprocessing_all_conjugated_systems!(mol::AbstractMolecule)
    mol_graph = mol.properties["mol_graph"]
    group_atomtypes_dict = Dict(1 => ["cc", "ce", "cp", "nc", "ne", "pc", "pe"],
                                2 => ["cd", "cf", "cq", "nd", "nf", "pd", "pf"])

    group1_conjugated_atoms = mol.atoms[(in(group_atomtypes_dict[1]).(mol.atoms[!,:atomtype])),:]

    group1_atoms_with_1_neighbor_in_group1_atoms = filter(x -> 
                countmap(in(group1_conjugated_atoms.number).(neighbors(mol_graph, x.number)))[true] == 1, 
                group1_conjugated_atoms)

    seen_atoms_vec = Vector{Int}()
    
    for conjAtm in eachrow(group1_atoms_with_1_neighbor_in_group1_atoms)
        if conjAtm.number in seen_atoms_vec
            continue
        end
        neighbors_in_conjugated_atoms = mol.atoms[(in(group_atomtypes_dict[1]).(mol.atoms[!,:atomtype]) .&& 
                                            in(neighbors(mol_graph, conjAtm.number)).(mol.atoms[!,:number])),:]
        
        depth = 1
        next_neigh_vec = neighbors_in_conjugated_atoms.number
        start_atom = conjAtm.number
        all_paths_from_curr_conjAtm = conjAtm_path_builder!(group1_conjugated_atoms, [start_atom, next_neigh_vec[1]], mol, depth)
        append!(seen_atoms_vec, reduce(vcat, all_paths_from_curr_conjAtm))
    end
end


function conjAtm_path_builder!(group1_conjugated_atoms::DataFrame, path_vec::Vector{Int}, mol::AbstractMolecule, depth::Int)
    mol_graph = mol.properties["mol_graph"]
    group_atomtypes_dict = Dict(1 => ["cc", "ce", "cp", "nc", "ne", "pc", "pe"],
                                2 => ["cd", "cf", "cq", "nd", "nf", "pd", "pf"])
    source_atom = path_vec[1]
    previous_atom = path_vec[lastindex(path_vec)-1]
    curr_atom = path_vec[lastindex(path_vec)]

    combination_df = DataFrame(:a1 => [previous_atom, curr_atom], :a2 => [curr_atom, previous_atom])
    bond_in_conjugated_atoms = innerjoin(combination_df, mol.bonds, on = [:a1, :a2])
    
    if bond_in_conjugated_atoms.order[1] == BondOrder.T(1)
        if mol.atoms.atomtype[previous_atom] in group_atomtypes_dict[2]
            mol.atoms.atomtype[curr_atom] = group_atomtypes_dict[2][findfirst(x -> x == mol.atoms.atomtype[curr_atom], group_atomtypes_dict[1])]
        end
    elseif bond_in_conjugated_atoms.order[1] == BondOrder.T(2)
        if mol.atoms.atomtype[previous_atom] in group_atomtypes_dict[1]
            mol.atoms.atomtype[curr_atom] = group_atomtypes_dict[2][findfirst(x -> x == mol.atoms.atomtype[curr_atom], group_atomtypes_dict[1])]
        end
    end
    next_neighbors_vec = filter(x -> !in(path_vec).(x) && x in neighbors(mol_graph, curr_atom) && 
                                    in(group1_conjugated_atoms.number).(x), neighborhood(mol_graph, source_atom, depth+1))
                                    
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