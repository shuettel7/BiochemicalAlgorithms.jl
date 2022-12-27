
export gaff_postprocessing_all_conjugated_systems!

function gaff_postprocessing_all_conjugated_systems!(mol::AbstractMolecule)
    mol_graph = mol.properties["mol_graph"]
    group_atomtypes_dict = Dict(1 => ["cc", "ce", "cp", "nc", "ne", "pc", "pe"],
                                2 => ["cd", "cf", "cq", "nd", "nf", "pd", "pf"])

    conjugated_atoms = mol.atoms[(in(group_atomtypes_dict[1]).(mol.atoms[!,:atomtype])),:]
    altered_atoms_vec = Vector{Int}()
    seen_atoms_vec = Vector{Int}()
    for conjAtm in eachrow(conjugated_atoms)
        neighbors_in_conjugated_atoms = mol.atoms[(in(group_atomtypes_dict[1]).(mol.atoms[!,:atomtype]) .&& 
                                            in(neighbors(mol_graph, conjAtm.number)).(mol.atoms[!,:number])),:]
        combination_df = DataFrame(:a1 => vcat(neighbors_in_conjugated_atoms.number, repeat([conjAtm.number],nrow(neighbors_in_conjugated_atoms))),
                            :a2 => vcat(repeat([conjAtm.number], nrow(neighbors_in_conjugated_atoms)), neighbors_in_conjugated_atoms.number))
        bonds_in_conjugated_atoms = innerjoin(combination_df, mol.bonds, on = [:a1, :a2])
        
        if nrow(neighbors_in_conjugated_atoms) == 1
            mol.atoms.atomtype[neighbors_in_conjugated_atoms.number[1]] = bonds_in_conjugated_atoms.order[1] == BondOrder.Double ? 
                                        group_atomtypes_dict[2][findfirst(x -> x == neighbors_in_conjugated_atoms.atomtype[1], group_atomtypes_dict[1])] :
                                        neighbors_in_conjugated_atoms.atomtype[1]
        elseif nrow(neighbors_in_conjugated_atoms) == 2
            conjugated_double_bonded_neighbor = conjugated_double_bond.a1[1] == conjAtm.number[1] ? conjugated_double_bond.a2[1] : conjugated_double_bond.a1[1]
            mol.atoms.atomtype[conjAtm.number[1]] = bonds_in_conjugated_atoms ### hier
                            group_atomtypes_dict[2][findfirst(x -> x == mol.atoms.atomtype[conjugated_double_bonded_neighbor], group_atomtypes_dict[1])]
        end
        println(mol.atoms.atomtype)
        readline()
        println(neighbors_in_conjugated_atoms)
    end
end