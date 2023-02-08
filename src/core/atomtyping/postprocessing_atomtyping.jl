
export gaff_postprocessing_all_conjugated_systems!


function gaff_postprocessing_all_conjugated_systems!(mol::Molecule)
    mol_graph = mol.properties["mol_graph"]
    group_atomtypes_dict = Dict{Int, Vector{String}}(
                                1 => ["cc", "ce", "cg", "cp", "nc", "ne", "pc", "pe"],
                                2 => ["cd", "cf", "ch", "cq", "nd", "nf", "pd", "pf"])

    group1_conjugated_atoms = mol.atoms[(in(group_atomtypes_dict[1]).(mol.atoms[!,:atomtype])),:]

    seen_atoms_vec = Vector{Int}()
    seen_neigh_atoms_vec = Vector{Int}()
    
    for conjAtm in eachrow(group1_conjugated_atoms)
        append!(seen_atoms_vec, [conjAtm.number[1]])

        # adapt the current atom to the neighboring atoms that have already been changed
        seen_next_neigh_conjugated_atoms = mol.atoms[((in(seen_atoms_vec).(mol.atoms[!,:number]) .|| 
                                    in(seen_neigh_atoms_vec).(mol.atoms[!,:number])) .&&
                                    (in(group_atomtypes_dict[1]).(mol.atoms[!,:atomtype]) .||
                                    in(group_atomtypes_dict[2]).(mol.atoms[!,:atomtype])) .&& 
                                    in(neighbors(mol_graph, conjAtm.number)).(mol.atoms[!,:number])),:]

        seen_next_neigh_vec = sort(seen_next_neigh_conjugated_atoms.number)
        
        if lastindex(seen_next_neigh_vec) > 0
            previous_atom = conjAtm.number[1]
            curr_atom = seen_next_neigh_vec[1]
            combination_df = DataFrame(:a1 => [previous_atom, curr_atom], :a2 => [curr_atom, previous_atom])
            bond_in_conjugated_atoms = innerjoin(combination_df, mol.bonds, on = [:a1, :a2])

            if haskey(bond_in_conjugated_atoms.properties[1], "TRIPOS_tag") && bond_in_conjugated_atoms.properties[1]["TRIPOS_tag"] == "ar"
                if mol.atoms.atomtype[curr_atom] in group_atomtypes_dict[1] && in(mol.atoms.atomtype[previous_atom], group_atomtypes_dict[1])
                    mol.atoms.atomtype[previous_atom] = group_atomtypes_dict[2][findfirst(x -> x == mol.atoms.atomtype[previous_atom], group_atomtypes_dict[1])]
                end
            elseif bond_in_conjugated_atoms.order[1] == BondOrder.T(1)
                if mol.atoms.atomtype[curr_atom] in group_atomtypes_dict[2] && in(mol.atoms.atomtype[previous_atom], group_atomtypes_dict[1])
                    mol.atoms.atomtype[previous_atom] = group_atomtypes_dict[2][findfirst(x -> x == mol.atoms.atomtype[previous_atom], group_atomtypes_dict[1])]
                end
            elseif bond_in_conjugated_atoms.order[1] == BondOrder.T(2)
                if mol.atoms.atomtype[curr_atom] in group_atomtypes_dict[1] && in(mol.atoms.atomtype[previous_atom], group_atomtypes_dict[1])
                    mol.atoms.atomtype[previous_atom] = group_atomtypes_dict[2][findfirst(x -> x == mol.atoms.atomtype[previous_atom], group_atomtypes_dict[1])]
                end
            elseif bond_in_conjugated_atoms.order[1] == BondOrder.T(3)
                if mol.atoms.atomtype[curr_atom] in group_atomtypes_dict[1] && in(mol.atoms.atomtype[previous_atom], group_atomtypes_dict[1])
                    mol.atoms.atomtype[previous_atom] = group_atomtypes_dict[2][findfirst(x -> x == mol.atoms.atomtype[previous_atom], group_atomtypes_dict[1])]
                end
            end
        end

        # now adapt the previously unseen neighbors to the current atom
        neighbors_in_conjugated_atoms = mol.atoms[((in(group_atomtypes_dict[1]).(mol.atoms[!,:atomtype]) .||
                                            in(group_atomtypes_dict[2]).(mol.atoms[!,:atomtype])) .&& 
                                            (.!in(seen_atoms_vec).(mol.atoms[!,:number]) .&&
                                                .!in(seen_neigh_atoms_vec).(mol.atoms[!,:number])) .&&
                                            in(neighbors(mol_graph, conjAtm.number)).(mol.atoms[!,:number])),:]
        
        next_neigh_vec = neighbors_in_conjugated_atoms.number

        for neigh in next_neigh_vec
            previous_atom = conjAtm.number[1]
            curr_atom = neigh
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
            elseif bond_in_conjugated_atoms.order[1] == BondOrder.T(3)
                if mol.atoms.atomtype[previous_atom] in group_atomtypes_dict[1] && in(mol.atoms.atomtype[curr_atom], group_atomtypes_dict[1])
                    mol.atoms.atomtype[curr_atom] = group_atomtypes_dict[2][findfirst(x -> x == mol.atoms.atomtype[curr_atom], group_atomtypes_dict[1])]
                end
            end
            append!(seen_neigh_atoms_vec, [neigh])
        end
    end
end