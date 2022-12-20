using BiochemicalAlgorithms

export CES_parser, CES_processor


function CES_parser(colstring::String)
    
    open_chars_list = ['(', '[', '<']
    close_chars_list = [')', ']', '>']
    logic_chars = [','] # '.' not necessary because "or"-cases of CES-Atoms are in other rows of the DEF File. 
    # '.' logic is therefore only found in the APS which follows a CES-Atom and is processed with a function 
    # based on the APS_processor

    bracket_logic_list = findall(x -> x in vcat(open_chars_list,close_chars_list,logic_chars), colstring)
    
    layers_df = DataFrame([Vector{Int}(), Vector{Int}(), Vector{Union{Int,String}}(), 
                            Vector{String}(), Vector{Int}(), Vector{String}(), 
                            Vector{String}(), Vector{Int}(), Vector{Vector{Int}}()],
                        ["LayerId", "LayerDepth", "GenericName", "Element", "NumNeighbors", "ElementWithNeighborCount", 
                            "CES_APS", "SourceRowNum", "ContainsLayerId"])
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
            
            if layer > 1 
                owner_layer_df = layers_df[(layers_df.LayerDepth .== layer-1), :]
                layers_df.SourceRowNum[nrow(layers_df)] = owner_layer_df.LayerId[nrow(owner_layer_df)]
                push!(owner_layer_df.ContainsLayerId[nrow(owner_layer_df)], layers_df.LayerId[nrow(layers_df)])
            end
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
        push!(molecule_paths_vecs, path_checker(cesColdata, atmprops_df, mol, i, i, layerItem.LayerId, 1))
    end
    if !isempty(filter(x -> x > 1, values(countmap(molecule_paths_vecs)))) || !isempty(filter(x -> in(x).(false), molecule_paths_vecs))
        return false
    else
        return true
    end
end


function path_checker(cesColdata::DataFrame, atmprops_df::DataFrameRow, mol::AbstractMolecule, absSourceAtm::Int, relSourceAtm::Int, cesAtomId::Int, depth::Int)
    mol_graph = mol.properties["mol_graph"]

    # List of current neighbor atoms at certain depth/distance from source 
    curr_atm_neighbors = filter(x -> !(x in neighborhood(mol_graph, absSourceAtm, depth-1)) && x in neighbors(mol_graph, relSourceAtm), neighborhood(mol_graph, absSourceAtm, depth))
    
    # DataFrame with ID/number, Element, and number of neighbors of the currently spectated atoms
    atm_df = DataFrame("AtmNum" => relSourceAtm, "Element" => enumToString(mol.atoms.element[relSourceAtm]), "NumNeighbors" => lastindex(neighbors(mol_graph, relSourceAtm)))
    cesRow = cesColdata[cesAtomId,:]
    next_cesAtoms = cesRow.ContainsLayerId

    # Path, list of atoms in chain in molecule resembling the CES 
    path_lists = Vector{Vector{Union{Int, Bool}}}()

    # compare and filter by elements and number of neighbors
    if isempty(cesRow.Element) && !isempty(cesRow.NumNeighbors)
        atm_df = atm_df[(atm_df.Element .== cesRow.Element),:]
    elseif !isempty(cesRow.Element) && isempty(cesRow.NumNeighbors)
        atm_df = atm_df[(atm_df.NumNeighbors .== cesRow.NumNeighbors),:]
    elseif !isempty(cesRow.Element) && !isempty(cesRow.NumNeighbors)
        atm_df = atm_df[(atm_df.NumNeighbors .== cesRow.NumNeighbors .&& atm_df.Element .== cesRow.Element),:]
    end

    # compare and check the CES_APS properties against the potential atoms of the molecule
    if !isempty(atm_df)
        for atm_num in atm_df.AtmNum
            push!(path_lists, [])
            # if there are matching atoms, check if there is a demanded CES_APS
            CES_APS_check = CES_APS_processor(cesRow.CES_APS, atmprops_df, atm_num, relSourceAtm, mol)
            
            if CES_APS_check && isempty(next_cesAtoms)
                return [atm_num]
            elseif CES_APS_check && !isempty(next_cesAtoms)
                for (cesIter, nextCesAtom) in enumerate(next_cesAtoms)
                    # for each matching atom to a given ces_aps, add the atom.number to the path_list
                    if CES_APS_check && !isempty(next_cesAtoms)
                        append!(path_lists, path_checker(cesColdata, atmprops_df, mol, absSourceAtm, atm_num, nextCesAtom, depth+1))
                    elseif !CES_APS_check
                        append!(path_lists, false)
                    end
                end
            elseif !CES_APS_check
                return [false]
            end
        end
    end
    println("ping1")
    return filter(x -> !in(x).(false), path_lists)
    
    ### for BFS use this for loop here
    # for listnum in 1:lastindex(ces_potential_partners_list) 
    #     for nextAtm in ces_potential_partners_list[listnum]
    #         return path_checker(cesColdata, atmprops_df, mol, absSourceAtm, nextAtm, layerrow.ContainsLayerId, depth+1)
    #     end
    # end
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