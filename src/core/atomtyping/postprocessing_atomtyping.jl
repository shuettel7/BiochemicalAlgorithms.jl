
export gaff_postprocessing_all_conjugated_systems!

function gaff_postprocessing_all_conjugated_systems!(mol::AbstractMolecule)
    mol_graph = mol.properties["mol_graph"]
    group_atomtypes_dict = Dict(1 => ["cc", "ce", "cp", "nc", "ne", "pc", "pe"],
                                2 => ["cd", "cf", "cq", "nd", "nf", "pd", "pf"])

    cat_chem_cycles_array = reduce(vcat, mol.properties["chem_cycle_array"])

    group1_conjugated_atoms = mol.atoms[(in(group_atomtypes_dict[1]).(mol.atoms[!,:atomtype])),:]

    seen_atoms_vec = Vector{Int}()
    
    for conjAtm in eachrow(group1_conjugated_atoms)
        if conjAtm.number in seen_atoms_vec
            continue
        end
        neighbors_in_conjugated_atoms = mol.atoms[(in(group_atomtypes_dict[1]).(mol.atoms[!,:atomtype]) .&& 
                                            in(neighbors(mol_graph, conjAtm.number)).(mol.atoms[!,:number])),:]
        
        if nrow(neighbors_in_conjugated_atoms) == 1 && !in(seen_atoms_vec).(conjAtm.number[1])
            depth = 1
            next_neigh = neighbors_in_conjugated_atoms.number[1]
            previous_atom = conjAtm.number
            push!(seen_atoms_vec, conjAtm.number[1])
            while lastindex(next_neigh) > 0 
                combination_df = DataFrame(:a1 => vcat(next_neigh[1], previous_atom[1]), :a2 => vcat(previous_atom[1], next_neigh[1]))
                bond_in_conjugated_atoms = innerjoin(combination_df, mol.bonds, on = [:a1, :a2])
                if bond_in_conjugated_atoms.order[1] == BondOrder.T(1)
                    if mol.atoms.atomtype[previous_atom[1]] in group_atomtypes_dict[2]
                        mol.atoms.atomtype[next_neigh[1]] = group_atomtypes_dict[2][findfirst(x -> x == mol.atoms.atomtype[next_neigh[1]], group_atomtypes_dict[1])]
                    end
                elseif bond_in_conjugated_atoms.order[1] == BondOrder.T(2)
                    if mol.atoms.atomtype[previous_atom[1]] in group_atomtypes_dict[1]
                        mol.atoms.atomtype[next_neigh[1]] = group_atomtypes_dict[2][findfirst(x -> x == mol.atoms.atomtype[next_neigh[1]], group_atomtypes_dict[1])]
                    end
                end
                push!(seen_atoms_vec, next_neigh[1])
                println("next_neigh[1]:  ", next_neigh[1], "   previous_atom:  ", previous_atom[1])
                depth += 1
                previous_atom = copy(next_neigh)
                next_neigh = filter(x -> x in neighborhood(mol_graph, conjAtm.number, depth) &&
                                    !(x in neighborhood(mol_graph, conjAtm.number, depth-1)) && 
                                    !in(seen_atoms_vec).(x), group1_conjugated_atoms.number)
            end
        end
    end
end