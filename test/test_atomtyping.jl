
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


@testset "APS_processor" begin
    atmprops_df1 = DataFrameRow(DataFrame("BondTypes" => [["sb", "2sb", "db", "1db", "AR5", "RG3"]]), 1)
    atmprops_df2 = DataFrameRow(DataFrame("BondTypes" => [["SB", "2SB", "DB", "1DB", "NR"]]), 1)

    def_file_col_string1 = "[sb,AR5,RG3]"
    def_file_col_string2 = "[sb,db,AR5,RG3]"
    def_file_col_string3 = "[db,AR5,RG3.RG6.RG12]"
    def_file_col_string4 = "[DB,NR]"
    def_file_col_string5 = "[SB,DB,1DL,NR.RG3.RG6.RG12]"
    def_file_col_string6 = "[2SB,DB.DL.2DB]"

    @test APS_processor(def_file_col_string1, atmprops_df1) == true
    @test APS_processor(def_file_col_string2, atmprops_df1) == true    
    @test APS_processor(def_file_col_string3, atmprops_df1) == true
    @test APS_processor(def_file_col_string4, atmprops_df2) == true
    @test APS_processor(def_file_col_string5, atmprops_df2) == false
    @test APS_processor(def_file_col_string6, atmprops_df2) == true

end