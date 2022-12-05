using BiochemicalAlgorithms


function CES_parser(colstring::String, atmprops_df::DataFrameRow, mol::AbstractMolecule, layer::Int)
    if colstring[1] == '(' && colstring[lastindex(colstring)] == ')'
        colstring = colstring[2:lastindex(colstring)-1]
    end
    
    open_and_logic_chars_list = ['(', '[', '<', ',', '.']
    close_and_logic_chars_list = [')', ']', '>', ',', '.']

    open_round_brackets = count(==('('), colstring)
    open_edgy_brackets = count(==('['), colstring)
    commas_in_string = count(==(','), colstring)
    dots_in_string = count(==('.'), colstring)

    open_rounded_bracket_list = findall('(', colstring)
    close_rounded_bracket_list = findall(')', colstring)
    open_edgy_bracket_list = findall('[', colstring)
    close_edge_bracket_list = findall(']', colstring)
    comma_list = findall(',', colstring)
    dot_list = findall('.', colstring)

    next_open_rounded_bracket = findnext('(', colstring, 1)
    next_close_rounded_bracket = findnext(')', colstring, 1)
    next_open_edgy_bracket = findnext('[', colstring, 1)
    next_close_edge_bracket = findnext(']', colstring, 1)
    next_comma_list = findnext(',', colstring, 1)
    next_dot_list = findnext('.', colstring, 1)

    this_layer_expr = Expr(:||)
    next_layer_expr = Expr(:&&)
    
    sections_df = DataFrame([Vector{Int}(), Vector{String}(), Vector{Int}(),
                    Vector{Int}(), Vector{Vector{Int}}(), Vector{Vector{Int}}()],
                    ["NumOfType", "Type", "Open", "Closed", "Subgroups", "Intersections"])
    for (i,char) in enumerate(colstring)
        if char == '('
            layer += 1
            push!(next_layer_expr, CES_parser(colstring[i:lastindex(colstring)],atmprops_df,mol,layer))
        elseif char == ')'
            return 
        elseif char == '[' 
            ### Question about what eg. [sb'] means, Antechamber instructions unclear
            push!(sections_df, )
            return in(atmprops_df.BondTypes).()
        elseif char == ','


        elseif char == '.'

        end
    end 

    bool_check_expr = Expr(:&&)
    layer_list = Vector{Vector{current_CES_atom}}()

    while open_round_brackets > 0
        curr_CES_atom = colstring[1:min(next_rounded, next_edgy, next_comma)] 

    
end


function CES_processor()
    # use DEF file and information from CES_parser process to cycle through
    # DEF file until Atomtype is assigned or no more criteria fit (return error?)
end


mutable struct current_CES_atom
    element::String
    num_neigh::String
    ele_w_neigh_num::String
    atomic_property::Vector{String}()
    generic_name::String
end