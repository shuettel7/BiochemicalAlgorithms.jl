@testset "DEF_file_parser" begin
    def_file = load_atomtyping_DEF("data/TEST_DEF_file_parser_ATOMTYPE_GFF.DEF")

    @test def_file isa DataFrame
    @test nrow(def_file) == 296
    @test ncol(def_file) == 8
    @test def_file.type_name[9] == "ca"
    @test all(in(["*"]).(def_file.residue_names))
    @test typeof(def_file.atomic_number) == Vector{Int64}
    @test def_file.atomic_property[5] == "[1DB,0DL]"
    @test def_file.CES[18] == "(XD3[sb',db])"

end

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