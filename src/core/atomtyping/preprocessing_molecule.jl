
using BiochemicalAlgorithms
using Graphs, SimpleWeightedGraphs, StatsBase

export PreprocessingMolecule!, ClearPreprocessingMolecule!


function PreprocessingMolecule!(mol::AbstractMolecule)
    # Graph representation of Molecule
    mol_graph = build_graph(mol)
    adj_matrix = adjacency_matrix(mol_graph)
    mol_wgraph = build_weighted_graph(mol)
    wgraph_adj_matrix = adjacency_matrix(mol_wgraph)

    # Cycle detection and list
    chem_cycle_list = cycle_basis(mol_graph)
    ring_intersections_matrix = cycle_intersections(chem_cycle_list)
    ring_class_list = aromaticity_type_processor(chem_cycle_list, wgraph_adj_matrix, ring_intersections_matrix, mol)

    ### Assign properties to Atoms 
    ### TODO: Bonds, when properties are introduced to Bond NamedTuple (e.g. adding ar for aromatic, am for amide)
    for (j, sublist) in enumerate(chem_cycle_list)
        for k in sublist
            if !haskey(mol.atoms.properties[k], "CycleList")
                mol.atoms.properties[k]["CycleList"] = [sublist]
                mol.atoms.properties[k]["CycleSize"] = [string("RG", lastindex(sublist))]
            elseif !all(in(mol.atoms.properties[k]["CycleList"]).([sublist]))
                mol.atoms.properties[k]["CycleList"] = append!(mol.atoms.properties[k]["CycleList"], [sublist])
                mol.atoms.properties[k]["CycleSize"] = append!(mol.atoms.properties[k]["CycleSize"], [string("RG", lastindex(sublist))])
            end
        end
    end

    for i = (1:nrow(mol.atoms))
        mol.atoms.properties[i]["ElementWithNeighbourCount"] = string(enumToString(mol.atoms.element[i]), lastindex(neighbors(mol_graph, i)))
        mol.atoms.properties[i]["AromaticityType"] = ring_class_list[i]
        for l in values(countmap(adjacency_matrix(mol_wgraph)[i, neighbors(mol_wgraph, i)]))
            if !haskey(mol.atoms.properties[i], "BondTypes")
                mol.atoms.properties[i]["BondTypes"] = [enumToString(BondShortOrderType(l))]
            elseif !all(in(mol.atoms.properties[i]["BondTypes"]).([enumToString(BondShortOrderType(l))]))
                mol.atoms.properties[i]["BondTypes"] = append!(mol.atoms.properties[i]["BondTypes"], [enumToString(BondShortOrderType(l))])
            end
        end
    end
end


function ClearPreprocessingMolecule!(mol::AbstractMolecule)
    for i = (1:nrow(mol.atoms))
        if haskey(mol.atoms.properties[i], "CycleList")
            delete!(mol.atoms.properties[i], "CycleList")
            delete!(mol.atoms.properties[i], "CycleSize")
        end
        if haskey(mol.atoms.properties[i], "ElementWithNeighbourCount")
            delete!(mol.atoms.properties[i], "ElementWithNeighbourCount")
        end
        if haskey(mol.atoms.properties[i], "AromaticityType")
            delete!(mol.atoms.properties[i], "AromaticityType")
        end
        if haskey(mol.atoms.properties[i], "BondTypes")
            delete!(mol.atoms.properties[i], "BondTypes")
        end
    end
end


function aromaticity_type_processor(LList::Vector{Vector{Int64}}, wgraph_adj::Graphs.SparseMatrix, inters_matrix::Matrix{Vector{Int64}}, mol::AbstractMolecule)
    ring_class_list = Vector{Vector{String}}(undef, nrow(mol.atoms))
    for i = (1:nrow(mol.atoms))
        ring_class_list[i] = ["NG"]
    end
    for (numvlist, vlist) in enumerate(LList)

        # check if is O, N, or S present in Ring vertex list
        ONSP_present = false
        if true in in(mol.atoms.element[vlist]).([Elements.O,Elements.N,Elements.S,Elements.P])
            ONSP_present = true
        end
        
        # check number of pi electrons
        pi_elec = 0
        for bond = (1:lastindex(vlist)-1)
            if bond == 1 && Int(wgraph_adj[vlist[bond], vlist[lastindex(vlist)]]) == 2
                pi_elec += 2
            elseif Int(wgraph_adj[vlist[bond], vlist[bond+1]]) == 2
                pi_elec += 2
            end
        end
        if (pi_elec / lastindex(vlist)) == 1.0
            for x in vlist
                ring_class_list[x] = ["AR1"]
            end
        elseif (pi_elec / lastindex(vlist)) > 1/2 && (pi_elec / lastindex(vlist)) <= 1 && ONSP_present && lastindex(vlist) > 4
            for x in vlist
                ring_class_list[x] = ["AR2"]
            end
        elseif (pi_elec / lastindex(vlist)) > 1/2 && (pi_elec / lastindex(vlist)) < 1 && !ONSP_present && lastindex(vlist) > 4
            # check if Ring vlist has intersections with other rings in molecule and if these are aromatic
            for i = (1:lastindex(LList))
                if !isempty(inters_matrix[numvlist,i]) && lastindex(inters_matrix[numvlist, i]) == 2
                    atom1_bonds = filter(:a1 => n -> n == inters_matrix[numvlist,i][1], mol.bonds)
                    atom2_bonds = filter(:a1 => n -> n == inters_matrix[numvlist,i][2], mol.bonds)
                    if countmap(in([BondOrder.Double]).(atom1_bonds.order))[1] == 1 &&
                        countmap(in([BondOrder.Double]).(atom2_bonds.order))[1] == 1
                        for x in vlist
                            if !in(ring_class_list[x]).("AR1")
                                ring_class_list[x] = ["AR1"]
                            end
                        end 
                    end           
                end
            end
        elseif (pi_elec / lastindex(vlist)) < 1/2
            for x in vlist
                if !in(ring_class_list[x]).("AR5") && !in(ring_class_list[x]).("AR1")
                    ring_class_list[x] = ["AR5"]    
                end
            end
        end
    end
    return ring_class_list
end


function cycle_intersections(LList::Vector{Vector{Int64}})
    inters_matrix = Matrix{Vector{Int64}}(undef, lastindex(LList), lastindex(LList))
    for i = (1:lastindex(LList))
        curr1_list = copy(LList[i])
        for j = (1:lastindex(LList))
            curr2_list = copy(LList[j])
            if curr1_list != curr2_list
                inters_matrix[i,j] = inters_matrix[j,i] = intersect(curr1_list, curr2_list)
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