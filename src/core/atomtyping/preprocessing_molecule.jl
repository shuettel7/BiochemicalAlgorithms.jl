using BiochemicalAlgorithms
using Graphs, SimpleWeightedGraphs, StatsBase, EnumX, DataFramesMeta

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
        ring_atoms = DataFrame(:number => sublist)
        ring_bonds = DataFrame(:a1 => sublist, :a2 => append!(sublist[2:lastindex(sublist)], sublist[1]))
        if all(in(["AR1"]).(ring_class_list[sublist]))
            ring_bonds = DataFrame(:a1 => sublist, :a2 => append!(sublist[2:lastindex(sublist)], sublist[1]))
            @with innerjoin(ring_bonds, mol.bonds, on = [:a1, :a2]) begin
                push!.(:properties, "TRIPOS_tag" => "ar")
            end
            @with innerjoin(ring_bonds, mol.bonds, on = [:a1 => :a2, :a2 => :a1]) begin
                push!.(:properties, "TRIPOS_tag" => "ar")
            end
        end
        @with innerjoin(ring_atoms, mol.atoms, on = [:number]) begin
                push!.(:properties, "CycleListNum" => [j])
                push!.(:properties, "CycleSize" => [lastindex(sublist)])
        end
    end

    ElemWNeighbourCount_vector = @with mol.atoms @byrow string(enumToString(:element), lastindex(neighbors(mol_graph, :number)))
    @with mol.atoms @byrow push!(:properties, 
                        "ElementWithNeighborCount" => string(enumToString(:element), lastindex(neighbors(mol_graph, :number))),
                        "Neighbors" => neighbors(mol_graph, :number), 
                        "AromaticityType" => mol.properties["ring_class_list"][:number],
                        "SecondaryNeighbors" => secondary_neighbors(mol_graph, :number),
                        "BondTypes" => bondShortOrder_types(:number, mol_graph, wgraph_adj_matrix))        
    
    amide_list = amide_processor(mol, mol_graph, ElemWNeighbourCount_vector)
    if !isempty(amide_list)
        amide_bonded_atoms_df = DataFrame(NamedTuple{(:a1, :a2)}.(amide_list))
        @with innerjoin(amide_bonded_atoms_df, mol.bonds, on = [:a1, :a2]) begin
            push!.(:properties, "TRIPOS_tag" => "am")
        end
        @with innerjoin(amide_bonded_atoms_df, mol.bonds, on = [:a1 => :a2, :a2 => :a1]) begin
            push!.(:properties, "TRIPOS_tag" => "am")
        end
    end
end


function bondShortOrder_types(num::Int, mol_graph::Graph, wgraph_adj_matrix::Graphs.SparseMatrix)
    BondShortVec = Vector{String}()
    for i in keys(countmap(wgraph_adj_matrix[num, neighbors(mol_graph, num)]))
        push!(BondShortVec, enumToString(BondShortOrderType(Int(i))))
    end
    return BondShortVec
end


function secondary_neighbors(mol_graph::Graph, num::Int)
    SecNeighVec = Vector{Vector{Int}}()
    for (k, neigh) in enumerate(neighbors(mol_graph, num))
        append!(SecNeighVec, [[]])
        for sec_neigh in neighbors(mol_graph, neigh)
            if sec_neigh != num
                push!(SecNeighVec[k], sec_neigh)
            end
        end
    end
    return SecNeighVec
end


function amide_processor(mol::AbstractMolecule, mol_graph::Graph, ElementWNeighbourCount_vector::Vector{String})
    amide_bond_vector = Vector{Tuple{Int64, Int64}}()
    nitrogen_atoms_list = mol.atoms[(mol.atoms[!, :element] .== Elements.N), :number]
    for nitrogen in nitrogen_atoms_list
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
                        "weighted_graph_adj_matrix", "chem_cycle_list", 
                        "ring_intersections_matrix", "ring_class_list"]
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


function create_atom_preprocessing_df(mol::AbstractMolecule)
    # Create DataFrame for better accessibility and handling of atom properties 
    if !haskey(mol.atoms.properties[1], "ElementWithNeighborCount")
        println("Running PreprocessingMolecule! on $(mol.name)")
        ClearPreprocessingMolecule!(mol)
        PreprocessingMolecule!(mol)
    end

    col_names = ["AromaticityType", "ElementWithNeighborCount", "Neighbors", "SecondaryNeighbors", 
                 "BondTypes", "CycleSize", "CycleListNum"]
    atom_props_df = DataFrame([Vector{String}(), Vector{String}(), Vector{Vector{Int64}}(), 
                            Vector{Vector{Vector{Int64}}}(), Vector{Vector{String}}(), Vector{Vector{Int64}}(), 
                            Vector{Vector{Int64}}()], col_names)
    # mol_atom_props = DataFrame()

    # ### TODO: Metaprogramming version below?
    # @with innerjoin(DataFrame(mol.atoms.properties), atom_props_df) begin
    #     push.(mol_atom_props)
    # end
    # println(mol_atom_props)
    
    for (k, atm) in enumerate(eachrow(mol.atoms))
        push!(atom_props_df, (string(), string(), Vector{Int64}(), 
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
    ring_class_list = Vector{String}(undef, nrow(mol.atoms))
    for i = (1:nrow(mol.atoms))
        ring_class_list[i] = "NG"
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
                ring_class_list[x] = "AR1"
            end
        elseif (pi_elec / lastindex(vlist)) > 1/2 && (pi_elec / lastindex(vlist)) <= 1 && ONSP_present && lastindex(vlist) > 4
            for x in vlist
                ring_class_list[x] = "AR2"
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
                            if ring_class_list[x] != "AR1"
                                ring_class_list[x] = "AR1"
                            end
                        end 
                    end           
                end
            end
        elseif (pi_elec / lastindex(vlist)) < 1/2
            for x in vlist
                if ring_class_list[x] != "AR5" && ring_class_list[x] != "AR1"
                    ring_class_list[x] = "AR5"  
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


@enumx BondShortOrder begin
    sb = 1
    db = 2
    tb = 3 
    qb = 4
    un = 100
end

const BondShortOrderType = BondShortOrder.T
