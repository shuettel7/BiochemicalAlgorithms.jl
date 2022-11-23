
using BiochemicalAlgorithms
using Graphs, SimpleWeightedGraphs, StatsBase

export PreprocessingMolecule!, ClearPreprocessingMolecule!, create_atom_preprocessing_df


function PreprocessingMolecule!(mol::AbstractMolecule)
    # Graph representation of Molecule
    mol.properties["mol_graph"] = mol_graph = build_graph(mol)
    mol.properties["adjacency_matrix"] = adj_matrix = adjacency_matrix(mol.properties["mol_graph"])
    mol.properties["mol_weighted_graph"] = mol_wgraph = build_weighted_graph(mol)
    mol.properties["weighted_graph_adj_matrix"] = wgraph_adj_matrix = adjacency_matrix(mol.properties["mol_weighted_graph"])

    # Cycle detection and list
    mol.properties["chem_cycle_list"] = chem_cycle_list = cycle_basis(mol_graph)
    ring_intersections_matrix = Matrix{Vector{Int64}}(undef, lastindex(chem_cycle_list), lastindex(chem_cycle_list))
    if lastindex(chem_cycle_list) > 1
        mol.properties["ring_intersections_matrix"] = ring_intersections_matrix = cycle_intersections(chem_cycle_list)
    end
    mol.properties["ring_class_list"] = ring_class_list = aromaticity_type_processor(chem_cycle_list, wgraph_adj_matrix, ring_intersections_matrix, mol)

    ### Assign properties to Atoms and Bonds
    for (j, sublist) in enumerate(chem_cycle_list)
        for (k, sublitem) in enumerate(sublist)
            ### assign "ar" property in mol.bonds.properties["TRIPOS_tag] position of Dictionary
            if k != 1 
                if in(["AR1", "AR2", "AR3"]).(ring_class_list[sublitem][1]) &&
                    in(["AR1", "AR2", "AR3"]).(ring_class_list[sublist[k-1]][1])
                    bond_row = get_bond_row(mol, sublitem, sublist[k-1])
                    if bond_row != 0
                        mol.bonds.properties[bond_row]["TRIPOS_tag"] = "ar"
                    end
                end 
            else
                if in(["AR1", "AR2", "AR3"]).(ring_class_list[sublitem][1]) &&
                    in(["AR1", "AR2", "AR3"]).(ring_class_list[sublist[lastindex(sublist)]][1])
                    bond_row = get_bond_row(mol, sublitem, sublist[lastindex(sublist)])
                    if bond_row != 0
                        mol.bonds.properties[bond_row]["TRIPOS_tag"] = "ar"
                    end
                end
            end
            if !haskey(mol.atoms.properties[sublitem], "CycleListNum")
                mol.atoms.properties[sublitem]["CycleListNum"] = [j]
                mol.atoms.properties[sublitem]["CycleSize"] = [lastindex(sublist)]
            elseif !all(in(mol.atoms.properties[sublitem]["CycleListNum"]).([j]))
                mol.atoms.properties[sublitem]["CycleListNum"] = append!(mol.atoms.properties[sublitem]["CycleListNum"], [j])
                mol.atoms.properties[sublitem]["CycleSize"] = append!(mol.atoms.properties[sublitem]["CycleSize"], [lastindex(sublist)])
            end
        end
    end

    ElemWNeighbourCount_vector = Vector{String}()
    for i = (1:nrow(mol.atoms))
        mol.atoms.properties[i]["ElementWithNeighborCount"] = ElemWNeighbourCount = string(enumToString(mol.atoms.element[i]), lastindex(neighbors(mol_graph, i)))
        mol.atoms.properties[i]["Neighbors"] = AtomNeighbors = neighbors(mol_graph, i)
        mol.atoms.properties[i]["SecondaryNeighbors"] = Vector{Vector{Int}}()
        for (k, neigh) in enumerate(AtomNeighbors)
            append!(mol.atoms.properties[i]["SecondaryNeighbors"], [[]])
            for sec_neigh in neighbors(mol_graph, neigh)
                if sec_neigh != i
                    push!(mol.atoms.properties[i]["SecondaryNeighbors"][k], sec_neigh)
                end
            end
        end
        push!(ElemWNeighbourCount_vector, ElemWNeighbourCount)
        mol.atoms.properties[i]["AromaticityType"] = mol.properties["ring_class_list"][i]
        for l in values(countmap(wgraph_adj_matrix[i, neighbors(mol.properties["mol_weighted_graph"], i)]))
            if !haskey(mol.atoms.properties[i], "BondTypes")
                mol.atoms.properties[i]["BondTypes"] = [enumToString(BondShortOrderType(l))]
            elseif !all(in(mol.atoms.properties[i]["BondTypes"]).([enumToString(BondShortOrderType(l))]))
                mol.atoms.properties[i]["BondTypes"] = append!(mol.atoms.properties[i]["BondTypes"], [enumToString(BondShortOrderType(l))])
            end
        end
    end
    mol.properties["Amide_tag_list"] = amide_list = amide_processor(mol, mol_graph, wgraph_adj_matrix, ElemWNeighbourCount_vector)
    if !isempty(amide_list)
        for atompair in amide_list
            bond_row = get_bond_row(mol, atompair[1], atompair[2])
            if bond_row != 0 && !haskey(mol.bonds.properties[bond_row] ,"TRIPOS_tag")
                mol.bonds.properties[bond_row]["TRIPOS_tag"] = "am"
            end
        end
    else
        delete!(mol.properties, "Amide_tag_list")
    end
end


function amide_processor(mol::AbstractMolecule, mol_graph::SimpleGraph, wgraph_adj_matrix::Graphs.SparseMatrix, ElementWNeighbourCount_vector::Vector{String})
    amide_bond_vector = Vector{Tuple{Int64, Int64}}()
    nitrogen_atoms_list = mol.atoms[(mol.atoms[!, :element] .== Elements.N), :number]
    for i in nitrogen_atoms_list
        for j in neighbors(mol_graph, i)
            if lastindex(neighbors(mol_graph, j)) > 2 && in(ElementWNeighbourCount_vector[neighbors(mol_graph, j)]).("O1")
                push!(amide_bond_vector, (j,i))
            end
        end
    end
    return amide_bond_vector
end


function ClearPreprocessingMolecule!(mol::AbstractMolecule)
    if haskey(mol.properties, "mol_graph")
        delete!(mol.properties, "mol_graph")
    end
    if haskey(mol.properties, "adjacency_matrix")
        delete!(mol.properties, "adjacency_matrix")
    end
    if haskey(mol.properties, "mol_weighted_graph")
        delete!(mol.properties, "mol_weighted_graph")
    end
    if haskey(mol.properties, "weighted_graph_adj_matrix")
        delete!(mol.properties, "weighted_graph_adj_matrix")
    end
    if haskey(mol.properties, "chem_cycle_list")
        delete!(mol.properties, "chem_cycle_list")
    end
    if haskey(mol.properties, "ring_intersections_matrix")
        delete!(mol.properties, "ring_intersections_matrix")
    end
    if haskey(mol.properties, "ring_class_list")
        delete!(mol.properties, "ring_class_list")
    end
    if haskey(mol.properties, "Amide_tag_list")
        delete!(mol.properties, "Amide_tag_list")
    end
    for i = (1:nrow(mol.atoms))
        if haskey(mol.atoms.properties[i], "CycleListNum")
            delete!(mol.atoms.properties[i], "CycleListNum")
            delete!(mol.atoms.properties[i], "CycleSize")
        end
        if haskey(mol.atoms.properties[i], "ElementWithNeighborCount")
            delete!(mol.atoms.properties[i], "ElementWithNeighborCount")
        end
        if haskey(mol.atoms.properties[i], "AromaticityType")
            delete!(mol.atoms.properties[i], "AromaticityType")
        end
        if haskey(mol.atoms.properties[i], "BondTypes")
            delete!(mol.atoms.properties[i], "BondTypes")
        end
        if haskey(mol.atoms.properties[i], "Neighbors")
            delete!(mol.atoms.properties[i], "Neighbors")
        end
        if haskey(mol.atoms.properties[i], "SecondaryNeighbors")
            delete!(mol.atoms.properties[i], "SecondaryNeighbors")
        end
    end
    for i = (1:nrow(mol.bonds))
        if haskey(mol.bonds.properties[i], "TRIPOS_tag")
            delete!(mol.bonds.properties[i], "TRIPOS_tag")
        end
    end
end


function create_atom_preprocessing_df(mol::AbstractMolecule)
    # Create DataFrame for better accessibility and handling of atom properties 
    if !haskey(mol.atoms.properties[1], "ElementWithNeighborCount")
        println("Running PreprocessingMolecule! on $(mol.name)")
        ClearPreprocessingMolecule!(mol)
        PreprocessingMolecule!(mol)
    end

    col_names = ["AromaticityType", "ElementWithNeighborCount", "Neighbors", "SecondaryNeighbors", 
                 "BondTypes", "CycleSize", "CycleListNum"]
    atom_props_df = DataFrame([Vector{Vector{String}}(), Vector{String}(), Vector{Vector{Int64}}(), 
                            Vector{Vector{Vector{Int64}}}(), Vector{Vector{String}}(), Vector{Vector{Int64}}(), 
                            Vector{Vector{Int64}}()], col_names)
    for (k, atm) in enumerate(eachrow(mol.atoms))
        push!(atom_props_df, (Vector{String}(), "", Vector{Int64}(), 
                            Vector{Vector{Int64}}(), Vector{String}(), Vector{Int64}(), 
                            Vector{Int64}()))
        for col in col_names
            if haskey(atm.properties, col)
                atom_props_df[!,col][k] = atm.properties[col]
            end
        end
    end
    return atom_props_df
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