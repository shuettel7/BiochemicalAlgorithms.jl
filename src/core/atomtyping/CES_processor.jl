using BiochemicalAlgorithms


function CES_parser(colstring::String, atmprops_df::DataFrameRow, mol::AbstractMolecule, layer::Int, atom_num::Int)
    
    open_chars_list = ['(', '[', '<']
    close_chars_list = [')', ']', '>']
    logic_chars = [',', '.']

    bracket_logic_list = findall(x -> x in vcat(open_chars_list,close_chars_list,logic_chars), colstring)

    layer = 0
    layers_df = DataFrame([Vector{Int}(), Vector{String}(), Vector{Int}(), Vector{String}(), Vector{String}()],
                    ["LayerNum", "Element", "NumNeighbors", "ElementWithNeighborCount", "CES_APS" ])
    for (i,strindex) in enumerate(colstring[bracket_logic_list]) 
        substring = colstring[bracket_logic_list[i]+1:bracket_logic_list[i+1]-1]
        if strindex == '(' || (strindex == ',' && (colstring[bracket_logic_list[i]-1] == ')' || isnumeric(colstring[bracket_logic_list[i]-1])))
            if strinex == '('
                layer += 1
            end
            push!(layers_df, (layer, filter(!isnumeric, substring), filter(isnumeric, substring),string(Element, NumNeighbors), ""))
            if colstring[bracket_logic_list[i+1]] == '['
                layers_df.CES_APS[nrow(layers_df)] = colstring[bracket_logic_list[i+1]:findnext(']', colstring, bracket_logic_list[i+2])]
            end 
        elseif strindex == ')'
            layer -= 1
        end
    end

    layer_expr = Expr(:&&)
    
end


function CES_processor()
    # use DEF file and information from CES_parser process to cycle through
    # DEF file until Atomtype is assigned or no more criteria fit (return error?)
end


mutable struct CesAtom
    element::String
    num_neigh::String
    atomic_property::String
end