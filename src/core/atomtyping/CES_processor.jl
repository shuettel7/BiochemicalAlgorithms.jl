using BiochemicalAlgorithms

export CES_parser


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
                push!(owner_layer_df.ContainsRowNum[nrow(owner_layer_df)], layers_df.LayerId[nrow(layers_df)])
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


function CES_processor(cesColdata::DataFrame, atmprops_df::DataFrame, mol::AbstractMolecule, i::Int)
    # use DEF file and information from CES_parser process to cycle through
    # DEF file until Atomtype is assigned or no more criteria fit (return error?)

    XX_XA_XB_XD_dict = Dict("XX" => [Elements.C, Elements.N, Elements.O, Elements.S, Elements.P],
                            "XA" => [Elements.S, Elements.O],
                            "XB" => [Elements.N, Elements.P],
                            "XD" => [Elements.S, Elements.P])
    
    ### Path checker
    layer1_list = cesColdata[(cesColdata.LayerDepth .== layer[1]), :]
    depth = 1
    molecule_paths_vecs = path_checker(cesColdata, atmprops_df, mol, i, i, layer1_list, depth)
    
end


function path_checker(cesColdata::DataFrame, atmprops_df::DataFrame, mol::AbstractMolecule, absSourceAtm::Int, relSourceAtm::Int, cesAtomsList::Vector{Int}, depth::Int)
    mol_graph = mol.properties["mol_graph"]

    # List of potential next Atoms in chain
    ces_potential_partners_list = Vector{Vector{Int}}([[] for _ in 1:lastindex(cesAtomsList)])

    # List of current neighbor atoms at certain depth/distance from source 
    curr_atm_neighbors = filter(!(x -> x in neighborhood(mol_graph, absSourceAtm, depth-1)), neighborhood(mol_graph, absSourceAtm, depth))

    # DataFrame with ID/number, Element, and number of neighbors of the currently spectated atoms
    potential_atm_paths_df = DataFrame([Vector{Int}(), Vector{String}(), Vector{Union{Int,String}}()],["AtmNum","Element", "NumNeighbors"])
    for atm_num in curr_atm_neighbors
        push!(atom_neighbors_at_layer_depth_df, (atm_num, enumToString(mol.atoms.element[atm_num]), lastindex(neighbors(mol_graph, atm_num))))
    end

    # for each cesAtom in cesAtomsList, look for potential partners with fitting properties
    curr_cesColdata = copy(cesColdata[cesAtomsList,:])
    for (num,layerrow) in enumerate(eachrow(curr_cesColdata))
        curr_potential_atm_paths = copy(potential_atm_paths_df)

        # compare and filter by elements and number of neighbors
        if isempty(curr_cesColdata.Element) && !isempty(curr_cesColdata.NumNeighbors)
            curr_potential_atm_paths = curr_potential_atm_paths[(curr_potential_atm_paths.Element .== layerrow.Element),:]
        elseif !isempty(curr_cesColdata.Element) && isempty(curr_cesColdata.NumNeighbors)
            curr_potential_atm_paths = curr_potential_atm_paths[(curr_potential_atm_paths.NumNeighbors .== layerrow.NumNeighbors),:]
        elseif !isempty(curr_cesColdata.Element) && !isempty(curr_cesColdata.NumNeighbors)
            curr_potential_atm_paths = innerjoin(curr_potential_atm_paths, layerrow, on = [:Element, :NumNeighbors])
        end

        # compare and check the CES_APS properties against the potential atoms of the molecule
        CES_APS_bool = isempty(layerrow.CES_APS)
        if !empty(curr_potential_atm_paths)
            if !CES_APS_bool
                curr_potential_atm_nums = copy(curr_potential_atm_paths.AtmNum) 
                for atm_num in curr_potential_atm_nums
                    CES_APS_bool = CES_APS_processor(layerrow.CES_APS, atmprops_df[atm_num,:], atm_num, relSourceAtm, mol)
                    if CES_APS_bool
                        # for each matching atom to a given ces_aps, add the atom.number to the ces_potential_partners_list
                        push!(ces_potential_partners_list[num], atm_num)
                    end
                end
            else
                # if there are no CES_APS properties demanded, return the whole curr_potential_atm_paths.atm_num list
                ces_potential_partners_list[num] = curr_potential_atm_paths.atm_num
            end
        else
            # if there are no curr_potential_atm_paths, return false
            # this terminates the search, ultimately leading to no match, moving on to next def-File DataFrameRow 
            return false
        end
        ### for DFS use this for loop here
        for nextAtm in ces_potential_partners_list[num]
            return path_checker(cesColdata, atmprops_df, mol, absSourceAtm, nextAtm, layerrow.ContainsRowNum, depth+1)
        end
    end
    ### for BFS use this for loop here
    # for listnum in 1:lastindex(ces_potential_partners_list) 
    #     for nextAtm in ces_potential_partners_list[listnum]
    #         return path_checker(cesColdata, atmprops_df, mol, absSourceAtm, nextAtm, layerrow.ContainsRowNum, depth+1)
    #     end
    # end
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