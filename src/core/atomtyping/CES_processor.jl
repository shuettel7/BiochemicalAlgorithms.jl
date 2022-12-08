using BiochemicalAlgorithms

export CES_parser


function CES_parser(colstring::String)
    
    open_chars_list = ['(', '[', '<']
    close_chars_list = [')', ']', '>']
    logic_chars = [','] # '.' not necessary because "or"-cases of CES-Atoms are in other rows of the DEF File. 
    # '.' logic is therefore only found in the APS which follows a CES-Atom and is processed with a function 
    # based on the APS_processor

    bracket_logic_list = findall(x -> x in vcat(open_chars_list,close_chars_list,logic_chars), colstring)

    
    layers_df = DataFrame([Vector{Int}(), Vector{Int}(), Vector{String}(), Vector{Union{Int,String}}(), 
                            Vector{String}(), Vector{String}(), Vector{Int}(), Vector{Vector{Int}}()],
                        ["LayerId", "LayerDepth", "Element", "NumNeighbors", "ElementWithNeighborCount", 
                            "CES_APS", "SourceRowNum", "ContainsRowNum"])
    layer = 0
    layer_id = 0
    in_APS_bool = false
    for (i,strindex) in enumerate(colstring[bracket_logic_list]) 
        if strindex == '['
            in_APS_bool = true
        elseif strindex == ']' || strindex == '('
            in_APS_bool = false
        end
        if strindex == '(' || (strindex == ',' && !in_APS_bool)
            if strindex == '('
                layer += 1
            end
            layer_id += 1
            substring = colstring[bracket_logic_list[i]+1:bracket_logic_list[i+1]-1]
            curr_element = filter(!isnumeric, substring)
            if isempty(filter(isnumeric, substring)) || contains(substring, "dd")
                num_neighbors = ""
                if contains(substring, "dd")
                    curr_element = filter(!=('d'), substring)
                end
            else
                num_neighbors = parse(Int, filter(isnumeric, substring))
            end
            push!(layers_df, (layer_id, layer, curr_element , num_neighbors, string(curr_element, num_neighbors), "", 0, []))
            if colstring[bracket_logic_list[i+1]] == '['
                layers_df.CES_APS[nrow(layers_df)] = colstring[bracket_logic_list[i+1]:findnext(']', colstring, bracket_logic_list[i+2])]
            end 
            if layer > 1 
                owner_layer_df = layers_df[(layers_df.LayerDepth .== layer-1), :]
                layers_df.SourceRowNum[nrow(layers_df)] = owner_layer_df.LayerId[nrow(owner_layer_df)]
                push!(owner_layer_df.ContainsRowNum[nrow(owner_layer_df)], layers_df.LayerId[nrow(layers_df)])
            end
        elseif strindex == ')'
            layer -= 1
        end
    end
    return layers_df
end


function CES_processor(coldata::DataFrame, atmprops_df::DataFrame, mol::AbstractMolecule, i::Int)
    # use DEF file and information from CES_parser process to cycle through
    # DEF file until Atomtype is assigned or no more criteria fit (return error?)
end


function CES_APS_processor(colstring::String, atmprops_df::DataFrameRow, curr_atom::Int, pre_atom::Int, mol::AbstractMolecule)
    if colstring[1] == '[' && colstring[lastindex(colstring)] == ']'
        colstring = colstring[2:lastindex(colstring)-1]
    end
    and_list = split(colstring, ',')
    and_expr_conditionals = Expr(:&&)

    if !all(in(atmprops_df.BondTypes).(and_list))
        for andItem in and_list
            or_list = split(andItem, '.')
            or_expr_conditionals = Expr(:||)
            for orItem in or_list
                if contains(orItem, '\'')
                    push!(or_expr_conditionals.args, :(enumToString(BondShortOrderType(mol.properties["weighted_graph_adj_matrix"][curr_atom, pre_atom]) == orItem[1:lastindex(orItem)-1])))
                else
                    push!(or_expr_conditionals.args, :(in($(atmprops_df.BondTypes)).($orItem)))
                end
            end
            if eval(or_expr_conditionals) == false
                return false
            else
                push!(and_expr_conditionals.args, eval(or_expr_conditionals)) 
            end
        end
    else
        return true
    end
    return eval(and_expr_conditionals)
end