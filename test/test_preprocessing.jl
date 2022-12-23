
@testset "Preprocessing_Molecule" begin
    mol = load_pubchem_json("data/TEST_PREPROCESSING_MOLECULE_Efavirenz_Conformer3D_CID_64139.json")
    PreprocessingMolecule!(mol)
    @test length(mol.properties) == 9
    @test length(mol.atoms.properties) == count_atoms(mol)
    @test length(mol.bonds.properties) == count_bonds(mol)
    @test length(mol.atoms.properties[20]) == 7
    @test lastindex(mol.atoms.properties[13]["CycleListNum"]) == 2
    @test mol.atoms.properties[20]["Neighbors"] == [1, 17, 21]
    @test mol.atoms.properties[1]["SecondaryNeighbors"] == [17, 21]
    @test mol.bonds.properties[9]["TRIPOS_tag"] == "am"
    @test mol.bonds.properties[31]["TRIPOS_tag"] == "ar"
end

@testset "Clear_Preprocessing_Molecule" begin
    mol = load_pubchem_json("data/TEST_PREPROCESSING_MOLECULE_Efavirenz_Conformer3D_CID_64139.json")
    mol_props_names = ["mol_graph", "adjacency_matrix", "mol_weighted_graph", 
                        "weighted_graph_adj_matrix", "chem_cycle_list", 
                        "ring_intersections_matrix", "ring_class_list", "atmprops_df"]
    atom_props_names = ["CycleListNum", "CycleSize", "ElementWithNeighborCount",
                        "AromaticityType", "BondTypes", "Neighbors", "SecondaryNeighbors"]
    bond_props_names = ["TRIPOS_tag"]
    num_props_before_preprocessing = Dict("atoms" => length(mol.atoms.properties[20]), "mol" => length(mol.properties), "bonds" => length(mol.bonds.properties[1]))
    PreprocessingMolecule!(mol)
    ClearPreprocessingMolecule!(mol)
    
    @test !haskey(mol.properties, "ring_class_list")
    @test length(mol.atoms.properties[20]) == num_props_before_preprocessing["atoms"]
    for name in mol_props_names
        @test !haskey(mol.properties, name)
    end
    for num in 1:nrow(mol.atoms)
        for name in atom_props_names
            @test !haskey(mol.atoms.properties[num], name)
        end
    end
    for num in 1:nrow(mol.bonds)
        for name in bond_props_names
            @test !haskey(mol.bonds.properties[num], name)
        end
    end
end
