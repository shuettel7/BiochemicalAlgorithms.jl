@testset "Preprocessing_Molecule" begin
    mol = load_pubchem_json("data/TEST_PREPROCESSING_MOLECULE_Efavirenz_Conformer3D_CID_64139.json")
    PreprocessingMolecule!(mol)
    @test length(mol.properties) == 8
    @test length(mol.atoms.properties) == count_atoms(mol)
    @test length(mol.bonds.properties) == count_bonds(mol)
    @test length(mol.atoms.properties[20]) == 7
    @test mol.atoms.properties[20]["Neighbors"] == [1, 17, 21]
    @test mol.bonds.properties[9]["TRIPOS_tag"] == "am"
    @test mol.bonds.properties[31]["TRIPOS_tag"] == "ar"
end

@testset "Clear_Preprocessing_Molecule" begin
    mol = load_pubchem_json("data/TEST_PREPROCESSING_MOLECULE_Efavirenz_Conformer3D_CID_64139.json")
    num_props_before_preprocessing = length(mol.atoms.properties[20])
    PreprocessingMolecule!(mol)
    ClearPreprocessingMolecule!(mol)
    @test !haskey(mol.properties, "ring_class_list")
    @test length(mol.atoms.properties[20]) == num_props_before_preprocessing
    for i = (1:nrow(mol.atoms))
        @test !haskey(mol.atoms.properties[i], "CycleListNum")
    end
    for i = (1:nrow(mol.bonds))
        @test !haskey(mol.bonds.properties[i], "TRIPOS_tag")
    end
end
