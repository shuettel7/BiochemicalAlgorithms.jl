
using Graphs, SimpleWeightedGraphs, StatsBase, EnumX, DataFramesMeta

export PreprocessingMolecule!, ClearPreprocessingMolecule!


function PreprocessingMolecule!(mol::AbstractMolecule)
    ClearPreprocessingMolecule!(mol)
    # Graph representation of Molecule
    mol.properties["mol_graph"] = mol_graph = build_graph(mol)
    mol.properties["adjacency_matrix"] = adj_matrix = adjacency_matrix(mol.properties["mol_graph"])
    mol.properties["mol_weighted_graph"] = mol_wgraph = build_weighted_graph(mol)
    mol.properties["weighted_graph_adj_matrix"] = wgraph_adj_matrix = adjacency_matrix(mol.properties["mol_weighted_graph"])

    # Cycle detection and list
    mol.properties["chem_cycle_array"] = chem_cycle_array = cycle_basis(mol_graph)
    ring_intersections_matrix = Matrix{Vector{Int64}}(undef, lastindex(chem_cycle_array), lastindex(chem_cycle_array))
    if lastindex(chem_cycle_array) > 1
        mol.properties["ring_intersections_matrix"] = ring_intersections_matrix = cycle_intersections(chem_cycle_array)
    end
    mol.properties["atom_aromaticity_array"] = atom_aromaticity_array = atom_aromaticity_type_processor(chem_cycle_array, ring_intersections_matrix, mol)
    mol.properties["atom_conjugated_system_array"] = atom_conjugated_system_array = atom_conjugated_system_processor(chem_cycle_array, wgraph_adj_matrix, mol)

    ### Assign Cycle/Ring properties to Atoms and Bonds
    for (j, sublist) in enumerate(chem_cycle_array)
        ### add TRIPOS_tag to AR1 type bonds
        shifted_sublist = vcat(sublist[2:lastindex(sublist)], sublist[1])
        ring_bonds = DataFrame(:a1 => vcat(sublist, shifted_sublist), :a2 => vcat(shifted_sublist, sublist))
        if haskey(countmap(reduce(vcat, atom_aromaticity_array[sublist])), "AR1") && 
                countmap(reduce(vcat, atom_aromaticity_array[sublist]))["AR1"] == lastindex(sublist)
            @with innerjoin(ring_bonds, mol.bonds, on = [:a1, :a2]) begin
                push!.(:properties, "TRIPOS_tag" => "ar")
            end
        end
    
        ### add CycleListNumber and CycleSize to every ring atom property
        for num in sublist
            if !(haskey(mol.atoms.properties[num], "CycleListNum") && haskey(mol.atoms.properties[num], "CycleSize"))
                push!(mol.atoms.properties[num], "CycleListNum" => [j], "CycleSize" => [lastindex(sublist)])
            elseif !in(mol.atoms.properties[num]["CycleListNum"]).(j)
                append!(mol.atoms.properties[num]["CycleListNum"], [j])
                append!(mol.atoms.properties[num]["CycleSize"], [lastindex(sublist)])
            end
        end 
    end

    ### add further fields to atom properties
    ElemWNeighbourCount_vector = @with mol.atoms @byrow string(enumToString(:element), lastindex(neighbors(mol_graph, :number)))
    @with mol.atoms @byrow push!(:properties, 
                        "ElementWithNeighborCount" => string(enumToString(:element), lastindex(neighbors(mol_graph, :number))),
                        "Neighbors" => neighbors(mol_graph, :number), 
                        "AromaticityType" => mol.properties["atom_aromaticity_array"][:number],
                        "SecondaryNeighbors" => secondary_neighbors(mol_graph, :number),
                        "BondTypes" => bondShortOrder_types(:number, mol, mol_graph, wgraph_adj_matrix))        
    
    ### add amide TRIPOS_tag to amide bonds
    amide_array = amide_processor(mol, mol_graph, ElemWNeighbourCount_vector)
    if !isempty(amide_array)
        amide_bonded_atoms_df = DataFrame(NamedTuple{(:a1, :a2)}.(amide_array))
        append!(amide_bonded_atoms_df, NamedTuple{(:a2, :a1)}.(amide_array))
        @with innerjoin(amide_bonded_atoms_df, mol.bonds, on = [:a1, :a2]) begin
            push!.(:properties, "TRIPOS_tag" => "am")
        end
    end

    push!(mol.properties, "atmprops_df" => create_atom_preprocessing_df!(mol))
end


function atom_conjugated_system_processor(LList::Vector{Vector{Int64}}, wgraph_adj::Graphs.SparseMatrix, mol::AbstractMolecule)
    mol_graph = mol.properties["mol_graph"]
    all_cycle_atoms = Vector{Int}()
    conjugated_atoms_vec = isempty(LList) ? Vector{Int}() : reduce(vcat, LList)
    
    filtered_bonds_df = mol.bonds[(.!in(all_cycle_atoms).(mol.bonds.a1) .&& .!in(all_cycle_atoms).(mol.bonds.a2) .&& 
                                (mol.bonds.order .== BondOrder.Double .|| mol.bonds.order .== BondOrder.Triple)), :]
    possible_conjugated_atoms = keys(countmap(vcat(filtered_bonds_df.a1, filtered_bonds_df.a2)))
    for atmNum in possible_conjugated_atoms
        direct_bonds_countmap = countmap(mol.properties["weighted_graph_adj_matrix"][atmNum, neighbors(mol_graph, atmNum)])
        if in([Elements.C, Elements.N,Elements.O, Elements.S, Elements.P]).(mol.atoms.element[atmNum]) && 
            ((haskey(direct_bonds_countmap, 2.0) && direct_bonds_countmap[2.0] >= 1) ||
            (haskey(direct_bonds_countmap, 3.0) && direct_bonds_countmap[3.0] >= 1))
            number_of_conjugatable_neighbors = 0
            for neigh in neighbors(mol_graph, atmNum)
                neighbors_bonds_countmap = countmap(mol.properties["weighted_graph_adj_matrix"][neigh, neighbors(mol_graph, neigh)])
                if in([Elements.C, Elements.N, Elements.O, Elements.S, Elements.P]).(mol.atoms.element[neigh]) && 
                    ((haskey(neighbors_bonds_countmap, 2.0) && neighbors_bonds_countmap[2.0] >= 1) ||
                    (haskey(neighbors_bonds_countmap, 3.0) && neighbors_bonds_countmap[3.0] >= 1))
                    number_of_conjugatable_neighbors += 1
                end
            end
            if number_of_conjugatable_neighbors >= 2
                push!(conjugated_atoms_vec, atmNum)
            end
        end
    end
    return conjugated_atoms_vec
end


function bondShortOrder_types(num::Int, mol::AbstractMolecule, mol_graph::Graph, wgraph_adj_matrix::Graphs.SparseMatrix)
    BondShortVec = Vector{String}()
    bonds_dict = countmap(wgraph_adj_matrix[num, neighbors(mol_graph, num)])
    for i in keys(bonds_dict)
        curr_bond_str = enumToString(BondShortOrderType(Int(i)))
        if !in(mol.properties["atom_aromaticity_array"][num]).("NR") && !in(mol.properties["atom_aromaticity_array"][num]).("AR5") || in(mol.properties["atom_conjugated_system_array"]).(num)
            push!(BondShortVec, curr_bond_str, string(bonds_dict[i], curr_bond_str))
        else
            push!(BondShortVec, uppercase(curr_bond_str), string(bonds_dict[i], uppercase(curr_bond_str)))
        end
    end
    for prop in mol.properties["atom_aromaticity_array"][num]
        push!(BondShortVec, prop)
    end
    if in(mol.properties["atom_conjugated_system_array"]).(num)
        push!(BondShortVec, "DL")
    end
    return BondShortVec
end


function secondary_neighbors(mol_graph::Graph, num::Int)
    return filter!(!(x -> x in neighborhood(mol_graph, num, 1)), neighborhood(mol_graph, num, 2))
end


function amide_processor(mol::AbstractMolecule, mol_graph::Graph, ElementWNeighbourCount_vector::Vector{String})
    amide_bond_vector = Vector{Tuple{Int64, Int64}}()
    nitrogen_atoms_array = mol.atoms[(mol.atoms[!, :element] .== Elements.N), :number]
    for nitrogen in nitrogen_atoms_array
        for amide_neigh in neighbors(mol_graph, nitrogen)
            if lastindex(neighbors(mol_graph, amide_neigh)) > 2 && 
                in(ElementWNeighbourCount_vector[neighbors(mol_graph, amide_neigh)]).("O1")
                push!(amide_bond_vector, (amide_neigh,nitrogen))
            end
        end
    end
    return amide_bond_vector
end


function ClearPreprocessingMolecule!(mol::AbstractMolecule)
    mol_props_names = ["mol_graph", "adjacency_matrix", "mol_weighted_graph", 
                        "weighted_graph_adj_matrix", "chem_cycle_array", 
                        "ring_intersections_matrix", "atom_aromaticity_array", 
                        "atmprops_df", "atom_conjugated_system_array"]
    atom_props_names = ["CycleListNum", "CycleSize", "ElementWithNeighborCount",
                        "AromaticityType", "BondTypes", "Neighbors", "SecondaryNeighbors"]
    bond_props_names = ["TRIPOS_tag"]
    for name in mol_props_names
        delete!(mol.properties, name)
    end
    @with mol.atoms begin
        for name in atom_props_names
            delete!.(:properties, name)
        end
    end
    @with mol.bonds begin
        for name in bond_props_names
            delete!.(:properties, name)
        end
    end
end


function create_atom_preprocessing_df!(mol::AbstractMolecule)
    # Create DataFrame for better accessibility and handling of atom properties 
    
    col_names = ["AromaticityType", "ElementWithNeighborCount", "Neighbors", "SecondaryNeighbors", 
                 "BondTypes", "CycleSize", "CycleListNum"]
    atom_props_df = DataFrame([Vector{Vector{String}}(), Vector{String}(), Vector{Vector{Int64}}(), 
                            Vector{Vector{Int64}}(), Vector{Vector{String}}(), Vector{Vector{Int64}}(), 
                            Vector{Vector{Int64}}()], col_names)        
    
    for (k, atm) in enumerate(eachrow(mol.atoms))
        push!(atom_props_df, (Vector{String}(), string(), Vector{Int64}(), 
                        Vector{Int64}(), Vector{String}(), Vector{Int64}(), 
                        Vector{Int64}()))
        for col in intersect(col_names, keys(atm.properties))
            atom_props_df[!,col][k] = atm.properties[col]
        end
    end
    return atom_props_df
end


function atom_aromaticity_type_processor(LList::Vector{Vector{Int64}}, inters_matrix::Matrix{Vector{Int64}}, mol::AbstractMolecule)
    atom_ring_class_array = Vector{Vector{String}}(undef, nrow(mol.atoms))
    for i = (1:nrow(mol.atoms))
        if i in reduce(vcat, isempty(LList) ? Vector{Int}() : LList)
            atom_ring_class_array[i] = []
        else
            atom_ring_class_array[i] = ["NR"]
        end
    end
    for (numvlist, vlist) in enumerate(LList)

        # check if is O, N, or S present in Ring vertex list
        ONSP_present = false
        if true in in(mol.atoms.element[vlist]).([Elements.O,Elements.N,Elements.S,Elements.P])
            ONSP_present = true
        end
        
        # check number of pi electrons
        sublist = copy(vlist)
        shifted_sublist = vcat(sublist[2:lastindex(sublist)], sublist[1])
        ring_bonds = DataFrame(:a1 => [sublist; shifted_sublist], :a2 => [shifted_sublist; sublist])
        innerjoined_df = innerjoin(ring_bonds, mol.bonds, on = [:a1, :a2])
        pi_electrons = haskey(countmap(innerjoined_df.order), BondOrder.T(2)) ? countmap(innerjoined_df.order)[BondOrder.T(2)] * 2 : 0
        if (pi_electrons / lastindex(vlist)) == 1.0
            for x in vlist
                if !in(atom_ring_class_array[x]).("AR1") && !in(atom_ring_class_array[x]).(string("RG6", lastindex(vlist)))
                    prepend!(atom_ring_class_array[x], ["AR1"], [string("RG", lastindex(vlist)), string(1, "RG", lastindex(vlist))])
                elseif in(atom_ring_class_array[x]).("AR1") && in(atom_ring_class_array[x]).(string("RG6", lastindex(vlist)))
                    count_ring_index = findfirst(x -> "AR1", atom_ring_class_array[x])+2
                    atom_ring_class_array[x][count_ring_index] = string(parse(Int,atom_ring_class_array[x][count_ring_index][1])+1, "RG", lastindex(vlist))
                end
            end
        elseif (pi_electrons / lastindex(vlist)) > 1/2 && (pi_electrons / lastindex(vlist)) <= 1 && ONSP_present && lastindex(vlist) > 4
            for x in vlist
                if !in(atom_ring_class_array[x]).("AR2") && !in(atom_ring_class_array[x]).(string("RG", lastindex(vlist)))
                    append!(atom_ring_class_array[x], ["AR2"], [string("RG", lastindex(vlist)), string(1, "RG", lastindex(vlist))])
                elseif in(atom_ring_class_array[x]).("AR2") && in(atom_ring_class_array[x]).(string("RG", lastindex(vlist)))
                    count_ring_index = findfirst(x -> "AR2", atom_ring_class_array[x])+2
                    atom_ring_class_array[x][count_ring_index] = string(parse(Int,atom_ring_class_array[x][count_ring_index][1])+1, "RG", lastindex(vlist))
                end
            end
        elseif (pi_electrons / lastindex(vlist)) > 1/2 && (pi_electrons / lastindex(vlist)) < 1 && !ONSP_present && lastindex(vlist) > 4
            # check if Ring vlist has intersections with other rings in molecule and if these are aromatic
            for i = (1:lastindex(LList))
                if !isempty(inters_matrix[numvlist,i]) && lastindex(inters_matrix[numvlist, i]) == 2
                    atom_bonds = mol.bonds[(mol.bonds.a1 .== inters_matrix[numvlist,i][1] .|| 
                                            mol.bonds.a2 .== inters_matrix[numvlist,i][2]),:]
                    if haskey(countmap(atom_bonds.order), BondOrder.T(2)) && countmap(atom_bonds.order)[BondOrder.T(2)] == 1
                        for x in vlist
                            if !in(atom_ring_class_array[x]).("AR1") && !in(atom_ring_class_array[x]).(string("RG6", lastindex(vlist)))
                                prepend!(atom_ring_class_array[x], ["AR1"], [string("RG", lastindex(vlist)), string(1, "RG", lastindex(vlist))])
                            elseif in(atom_ring_class_array[x]).("AR1") && in(atom_ring_class_array[x]).(string("RG6", lastindex(vlist)))
                                count_ring_index = findfirst(x -> "AR1", atom_ring_class_array[x])+2
                                atom_ring_class_array[x][count_ring_index] = string(parse(Int,atom_ring_class_array[x][count_ring_index][1])+1, "RG", lastindex(vlist))
                            end
                        end
                    end           
                end
            end
        elseif (pi_electrons / lastindex(vlist)) < 1/2
            for x in vlist
                if !in(atom_ring_class_array[x]).("AR5") && !in(atom_ring_class_array[x]).(string("RG6", lastindex(vlist)))
                    append!(atom_ring_class_array[x], ["AR5", string("RG", lastindex(vlist)), string(1, "RG", lastindex(vlist))]) 
                elseif in(atom_ring_class_array[x]).("AR5") && in(atom_ring_class_array[x]).(string("RG6", lastindex(vlist)))
                    ring_prop_with_count_index = findfirst(x -> string("RG", lastindex(vlist)), atom_ring_class_array[x])+1
                    atom_ring_class_array[x][ring_prop_with_count_index] = string(parse(Int,atom_ring_class_array[x][ring_prop_with_count_index][1])+1, "RG", lastindex(vlist))
                end
            end
        end
    end
    return atom_ring_class_array
end


function cycle_intersections(LList::Vector{Vector{Int64}})
    inters_matrix = Matrix{Vector{Int64}}(undef, lastindex(LList), lastindex(LList))
    for i = (1:lastindex(LList))
        curr1_array = copy(LList[i])
        for j = (1:lastindex(LList))
            curr2_array = copy(LList[j])
            if curr1_array != curr2_array
                inters_matrix[i,j] = inters_matrix[j,i] = intersect(curr1_array, curr2_array)
            else
                inters_matrix[i,j] = inters_matrix[j,i] = []
            end
        end
    end
    return inters_matrix
end


function build_weighted_graph(mol::AbstractMolecule)
    mol_weighted_graph = SimpleWeightedGraph(nrow(mol.atoms))
    for i = (1:nrow(mol.bonds))
        add_edge!(mol_weighted_graph, mol.bonds.a1[i], mol.bonds.a2[i], Int(mol.bonds.order[i]))
    end
    return mol_weighted_graph
end


function build_graph(mol::AbstractMolecule)
    mol_graph = SimpleGraph(nrow(mol.atoms))
    for i = (1:nrow(mol.bonds))
        add_edge!(mol_graph, mol.bonds.a1[i], mol.bonds.a2[i])
    end
    return mol_graph
end


function enumToString(AnyEnum::Enum)
    return String(Symbol(AnyEnum))
end


@enumx BondShortOrder begin
    sb = 1
    db = 2
    tb = 3 
    qb = 4
    un = 100
end

const BondShortOrderType = BondShortOrder.T
