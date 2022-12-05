using BiochemicalAlgorithms


function CES_parser(colstring::String, atmprops_df::DataFrameRow, mol::AbstractMolecule, layer::Int, atom_num::Int)

    current_layer_neighbors = get_current_layer(mol.properties["mol_graph"], atom_num, )
    
    open_chars_list = ['(', '[', '<']
    close_chars_list = [')', ']', '>']
    logic_chars = [',', '.']

    bracket_logic_list = findall(x -> x in vcat(open_chars_list,close_chars_list,logic_chars), colstring)

    layer = 0
    layers_df = DataFrame([Vector{Int}(), Vector{CesAtom}()],
                    ["LayerNum","CesAtom"])
    for (i,strindex) in enumerate(colstring[bracket_logic_list]) 
        substring = colstring[bracket_logic_list[i]+1:bracket_logic_list[i+1]-1]
        if strindex == '('
            layer += 1
            push!(layers_df, (layer, CesAtom("", "", "")))
            layers_df.CesAtom.element = filter(!isnumeric, substring)
            layers_df.CesAtom.num_neigh = filter(isnumeric, substring)
        elseif strindex == ')'
            layer -= 1
        elseif strindex == '['
            layers_df.CesAtom.atomic_property[layer] = colstring[bracket_logic_list[i]:findnext(']', colstring, bracket_logic_list[i+1])]
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