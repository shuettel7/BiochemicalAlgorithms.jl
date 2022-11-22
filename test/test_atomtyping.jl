@testset "Preprocessing_Molecule" begin
    mol = load_pubchem_json("data/TEST_PREPROCESSING_MOLECULE_Efavirenz_Conformer3D_CID_64139.json")
    PreprocessingMolecule!(mol)
    @test length(mol.properties) == 8
    @test length(mol.atoms.properties) == 30
    @test length(mol.bonds.properties) == 32
    @test length(mol.atoms.properties[20]) == 7
    @test mol.atoms.properties[20]["Neighbors"] == [1, 17, 21]
    @test mol.bonds.properties[9]["TRIPOS_tag"] == "am"
end