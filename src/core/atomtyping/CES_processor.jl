using BiochemicalAlgorithms


function CES_parser(colstring::String, atmprops_df::DataFrameRow, mol::AbstractMolecule, atom_num::Int)
    
    open_chars_list = ['(', '[', '<']
    close_chars_list = [')', ']', '>']
    logic_chars = [','] # '.' not necessary because "or"-cases of CES-Atoms are in other rows of the DEF File. 
    # '.' logic is therefore only found in the APS which follows a CES-Atom and is processed with a function 
    # based on the APS_processor

    bracket_logic_list = findall(x -> x in vcat(open_chars_list,close_chars_list,logic_chars), colstring)

    layer = 0
    layers_df = DataFrame([Vector{Int}(), Vector{String}(), Vector{Int}(), Vector{String}(), Vector{String}(), Vector{Vector{Int}}()],
                    ["LayerNum", "Element", "NumNeighbors", "ElementWithNeighborCount", "CES_APS", "ContainsRowNum"])
    for (i,strindex) in enumerate(colstring[bracket_logic_list]) 
        substring = colstring[bracket_logic_list[i]+1:bracket_logic_list[i+1]-1]
        if strindex == '(' || (strindex == ',' && (colstring[bracket_logic_list[i]-1] == ')' || isnumeric(colstring[bracket_logic_list[i]-1])))
            if strinex == '('
                layer += 1
            end
            push!(layers_df, (layer, filter(!isnumeric, substring), parse(Int, filter(isnumeric, substring)), string(Element, NumNeighbors), "", []))
            if colstring[bracket_logic_list[i+1]] == '['
                layers_df.CES_APS[nrow(layers_df)] = colstring[bracket_logic_list[i+1]:findnext(']', colstring, bracket_logic_list[i+2])]
            end 
            if layer > 1 
                owner_layer_df = layers_df[(layers_df.LayerNum .== layer-1), :]
                push!(owner_layer_df.ContainsRowNum[nrow(owner_layer_df)], nrow(layers_df))
            end
        elseif strindex == ')'
            layer -= 1
        end
    end

    for rownum in 1:nrow(layers_df)
        rowNumsChain = Vector{Int}()
        layer_depth = 0
        copy_layer_df = copy(layers_df)
        while !isempty(copy_layers_df.ContainsRowNum[])
    end

    layer_expr = Expr(:&&)

    # for layernum in keys(countmap(layers_df.LayerNum))
    #     curr_layers_df = layers_df[(layers_df.LayerNum .== layernum),:]
    #     curr_layer_atoms = filter(!(x -> x in neighborhood(mol1.properties["mol_graph"], atom_num, layernum-1)), neighborhood(mol1.properties["mol_graph"], atom_num, layernum))
    #     for atomnum in curr_layer_atoms
    #         push!(layer_expr.args, )
    # end
    

end


function CES_processor()
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
                    push!(or_expr_conditionals.args, :($(enumToString(BondShortOrderType(mol.properties["weighted_graph_adj_matrix"][curr_atom, pre_atom])) == $orItem[1:lastindex(orItem)-1])))
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