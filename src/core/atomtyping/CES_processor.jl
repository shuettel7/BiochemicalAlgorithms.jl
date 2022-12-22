using BiochemicalAlgorithms

export CES_parser, CES_processor, path_builder


function CES_parser(colstring::String, mol::AbstractMolecule, sourceAtomNum::Int)
    
    open_chars_list = ['(', '[', '<']
    close_chars_list = [')', ']', '>']
    logic_chars = [','] # '.' not necessary because "or"-cases of CES-Atoms are in other rows of the DEF File. 
    # '.' logic is therefore only found in the APS which follows a CES-Atom and is processed with a function 
    # based on the APS_processor

    bracket_logic_list = findall(x -> x in vcat(open_chars_list,close_chars_list,logic_chars), colstring)
    
    layers_df = DataFrame([Vector{Int}(), Vector{Int}(), Vector{String}(), Vector{String}(),
                            Vector{Union{Int,String}}(), Vector{String}(), Vector{String}(),  
                            Vector{Union{Int,String}}(), Vector{Vector{Int}}()],
                        ["AtomId", "LayerDepth", "GenericName", "Element", "NumNeighbors", "ElementWithNeighborCount", 
                            "CES_APS", "SourceRowNum", "ContainsAtomId"])

    layer = 0
    layer_id = 0
    in_APS_bool = false
    sourceMolGraph = mol.properties["mol_graph"]
    sourceElemWNeighCount = mol.properties["atmprops_df"].ElementWithNeighborCount[sourceAtomNum]
    push!(layers_df, (layer, layer_id, "", enumToString(mol.atoms.element[sourceAtomNum]), 
            lastindex(neighbors(sourceMolGraph, sourceAtomNum)), 
            sourceElemWNeighCount, "", "", []))
    
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
                if curr_element[lastindex(curr_element)-1:lastindex(curr_element)] == "dd"
                    curr_element = chop(curr_element, tail = 2)
                end
            else
                num_neighbors = parse(Int, filter(isnumeric, substring))
            end
            ces_APS_substring = ""
            generic_name = ""
            push!(layers_df, (layer_id, layer, generic_name, curr_element, num_neighbors, 
                                string(curr_element, num_neighbors), ces_APS_substring, 0, []))
            
            owner_layer_df = layers_df[(layers_df.LayerDepth .== layer-1), :]
            layers_df.SourceRowNum[nrow(layers_df)] = owner_layer_df.AtomId[nrow(owner_layer_df)]
            push!(owner_layer_df.ContainsAtomId[nrow(owner_layer_df)], layers_df.AtomId[nrow(layers_df)])
        elseif strindex == ')'
            layer -= 1
        elseif strindex == '['
            layers_df.CES_APS[nrow(layers_df)] = colstring[bracket_logic_list[i]:findnext(']', colstring, bracket_logic_list[i+1])]
        elseif strindex == '<'
            layers_df.GenericName[nrow(layers_df)] = colstring[bracket_logic_list[i]+1:findnext('>', colstring, bracket_logic_list[i+1])-1]
        end
    end
    return layers_df
end


function CES_processor(cesColdata::DataFrame, atmprops_df::DataFrameRow, mol::AbstractMolecule, i::Int)
    # use DEF file and information from CES_parser process to cycle through
    # DEF file until Atomtype is assigned or no more criteria fit (return error?)

    XX_XA_XB_XD_dict = Dict("XX" => [Elements.C, Elements.N, Elements.O, Elements.S, Elements.P],
                            "XA" => [Elements.S, Elements.O],
                            "XB" => [Elements.N, Elements.P],
                            "XD" => [Elements.S, Elements.P])
    
    ### Path checker
    startlayer_list = DataFrame(cesColdata[(cesColdata.LayerDepth .== 1), :])
    molecule_paths_vecs = Vector{Vector{Union{Int, Bool}}}()
    for (iter, layerItem) in enumerate(eachrow(startlayer_list))
        push!(molecule_paths_vecs, path_checker(cesColdata, atmprops_df, mol, i, i, layerItem.AtomId, 1))
    end
    if !isempty(filter(x -> x > 1, values(countmap(molecule_paths_vecs)))) || !isempty(filter(x -> in(x).(false), molecule_paths_vecs))
        return false
    else
        return true
    end
end


function path_checker(CES_df::DataFrame, atmprops_df::DataFrameRow, previousCesAtomIds::Vector{Int}, atmPath_vec::Vector{Int}, mol::AbstractMolecule, absSourceAtm::Int, depth::Int)
    nextCesAtomId_list = CES_df[(CES_df.AtomId .== previousCesAtomIds[lastindex(previousCesAtomIds)]), :ContainsAtomId][1]
    mol_graph = mol.properties["mol_graph"]
    
    filtered_atm_paths = Vector{Vector{Int}}()

    # List of current neighbor atoms at certain depth/distance from source 
    curr_atm_neighbors = filter(x -> !(x in neighborhood(mol_graph, absSourceAtm, depth-1)) && 
                                x in neighbors(mol_graph, atmPath_vec[lastindex(atmPath_vec)]), 
                                neighborhood(mol_graph, absSourceAtm, depth))
    
    cesRow = CES_df[previousCesAtomIds[lastindex(previousCesAtomIds)],:]

    # create an expression that evaluates by elements and number of neighbors if demanded
    check_expr = Expr(:&&)
    # Element properties check
    if in(keys(XX_XA_XB_XD_dict)).(cesRow.Element)
        push!(check_expr.args, in(XX_XA_XB_XD_dict[cesRow.Element]).(enumToString(mol.atoms.element[atmPath_vec[lastindex[atmPath_vec]]])))
    elseif !isempty(cesRow.Element) 
        push!(check_expr.args, cesRow.Element == enumToString(mol.atoms.element[atmPath_vec[lastindex[atmPath_vec]]]))
    end
    # check number of neigbors
    if !isempty(cesRow.NumNeighbors)
        push!(check_expr.args, cesRow.NumNeighbors == lastindex(neighbors(mol_graph, atmPath_vec[lastindex(atmPath_vec)])))
    end
    # check CES_APS against that of current atom
    push!(check_expr.args, CES_APS_processor(cesRow.CES_APS, atmprops_df, atmPath_vec[lastindex[atmPath_vec]], 
                                                atmPath_vec[lastindex[atmPath_vec]-1], mol))
    # evaluate the built expression to a boolean
    check_Bool = eval(check_expr)

    # return if last there are no following links in CES
    if isempty(cesRow.ContainsAtomId) && check_Bool
        push!(filtered_atm_paths, atmPath_vec)
        return filtered_atm_paths
    end

    # start next instance of function if path is correct so far
    if check_Bool
        for atm_neigh in curr_atm_neighbors
            for nextCesAtom in nextCesAtomId_list
                append!(filtered_atm_paths, path_checker(CES_df, atmprops_df, vcat(previousCesAtomIds, [nextCesAtom]), vcat(atmPath_vec, [atm_neigh]), mol, absSourceAtm, depth+1))
            end
        end
    end
    return filtered_atm_paths
end


function CES_APS_processor(colstring::String, atmprops_df::DataFrameRow, curr_atom::Int, pre_atom::Int, mol::AbstractMolecule)
    if colstring[1] == '[' && colstring[lastindex(colstring)] == ']'
        colstring = colstring[2:lastindex(colstring)-1]
    end
    # to shorten path checker. if the CES_APS is empty then there are no demands for the properties
    if isempty(colstring)
        return true
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