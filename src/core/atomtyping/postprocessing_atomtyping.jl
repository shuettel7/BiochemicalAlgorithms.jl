
function gaff_postprocessing_all_conjugated_systems!(mol::AbstractMolecule)
    group_atomtypes_dict = Dict(1 => [cc, ce, cp, nc, ne, pc, pe],
                                2 => [cd, cf, cq, nd, nf, pd, pf])

    possible_conjugated_atoms = mol.atoms[(in.(group_atomtypes_dict[1]).(mol.atoms.atomtype)),:]
end