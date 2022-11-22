@testset "PDB" begin
        bpti = load_pdb("data/bpti.pdb")

        @test bpti.name == "bpti.pdb"

        @test count_atoms(bpti) == 454
        @test count_bonds(bpti) == 0

        pdb_5pti = load_pdb("data/5PTI.pdb")

        @test pdb_5pti.name == "5PTI.pdb"

        @test count_atoms(pdb_5pti) == 1087

        @test count_bonds(pdb_5pti) == 0

end

@testset "PubChem" begin
        mol = load_pubchem_json("data/aspirin_pug.json")
        
        @test mol.name == "data/aspirin_pug.json"

        @test count_atoms(mol) == 21
        @test count_bonds(mol) == 21


        # used for testing bond orders / file manually manipulated
        mol2 = load_pubchem_json("data/aspirin_pug_bonds.json")

        @test count_atoms(mol2) == 21
        @test count_bonds(mol2) == 21
end

@testset "mol2_import" begin
        mol = load_mol2("data/Import_test_sustiva_modified.mol2")
        
        @test mol.name == "data/Import_test_sustiva_modified.mol2"

        @test count_atoms(mol) == 21
        @test mol.atoms.name[21] == "C14"
        @test mol.bonds.properties[2]["TRIPOS_tag"] == "am"
        @test count_bonds(mol) == 23
        @test mol.bonds.order[2] == BondOrder.Single
        @test mol.bonds.properties[2]["TRIPOS_tag"] == "am"
        @test mol.bonds.order[12] == BondOrder.Unknown
        @test mol.bonds.properties[12]["TRIPOS_tag"] == "ar"
end

@testset "mol2_export" begin
        export_mol2(load_pubchem_json("data/Export_test_molecule_Sustiva_Efavirenz_Conformer3D_CID_64139.json"), "data/")

        @test readlines("data/Export_test_molecule_Sustiva_Efavirenz_Conformer3D_CID_64139.mol2")[1] == "@<TRIPOS>MOLECULE"
        @test readlines("data/Export_test_molecule_Sustiva_Efavirenz_Conformer3D_CID_64139.mol2")[8] == "@<TRIPOS>ATOM"
        @test readlines("data/Export_test_molecule_Sustiva_Efavirenz_Conformer3D_CID_64139.mol2")[39] == "@<TRIPOS>BOND"

        rm("data/Export_test_molecule_Sustiva_Efavirenz_Conformer3D_CID_64139.mol2")

        export_mol2(load_pdb("data/Export_test_molecule_6dny.pdb"), "data/")

        @test readlines("data/Export_test_molecule_6dny.mol2")[1] == "@<TRIPOS>MOLECULE"
        @test readlines("data/Export_test_molecule_6dny.mol2")[8] == "@<TRIPOS>ATOM"
        @test readlines("data/Export_test_molecule_6dny.mol2")[74] == "@<TRIPOS>SUBSTRUCTURE"

        rm("data/Export_test_molecule_6dny.mol2")
end