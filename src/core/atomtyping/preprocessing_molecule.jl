
using Graphs, SimpleWeightedGraphs, StatsBase, EnumX, DataFramesMeta, BiochemicalAlgorithms

export PreprocessingMolecule!, ClearPreprocessingMolecule!


function PreprocessingMolecule!(mol::Molecule)
    # Graph representation of Molecule
    mol.properties["mol_graph"] = mol_graph = build_graph(mol)
    mol.properties["adjacency_matrix"] = adjacency_matrix(mol.properties["mol_graph"])
    mol.properties["mol_weighted_graph"] = build_weighted_graph(mol)
    mol.properties["weighted_graph_adj_matrix"] = wgraph_adj_matrix = adjacency_matrix(mol.properties["mol_weighted_graph"])

    # Cycle detection and Vectors
    # Antechamber only uses rings of size 3-9, anything larger is non-ring (NG in paper, NR here)
    mol.properties["chem_cycle_array"] = chem_cycle_array = filter(x -> lastindex(x) <= 9, cycle_basis(mol_graph))
    mol.properties["atom_aromaticity_array"] = atom_aromaticity_array = atom_aromaticity_type_processor(chem_cycle_array, mol)
    mol.properties["atom_conjugated_system_array"] = atom_conjugated_system_processor(chem_cycle_array, mol)

    ### Assign Cycle/Ring properties to Atoms and Bonds
    for (j, subvec) in enumerate(chem_cycle_array)
        ### add TRIPOS_tag to AR1 type bonds
        shifted_subvec = vcat(subvec[2:lastindex(subvec)], subvec[1])
        ring_bonds = DataFrame(:a1 => vcat(subvec, shifted_subvec), :a2 => vcat(shifted_subvec, subvec))
        atoms_dict = countmap(reduce(vcat, atom_aromaticity_array[subvec]), alg=:dict)
        if haskey(atoms_dict, "AR1") && atoms_dict["AR1"] == lastindex(subvec)
            @with innerjoin(ring_bonds, mol.bonds, on = [:a1, :a2]) begin
                push!.(:properties, "TRIPOS_tag" => "ar")
            end
        end
    
        ### add CycleVectorNumber and CycleSize to every ring atom property
        for num in subvec
            if !(haskey(mol.atoms.properties[num], "CycleVectorNum") && haskey(mol.atoms.properties[num], "CycleSize"))
                push!(mol.atoms.properties[num], "CycleVectorNum" => [j], "CycleSize" => [lastindex(subvec)])
            elseif !in(j, mol.atoms.properties[num]["CycleVectorNum"])
                append!(mol.atoms.properties[num]["CycleVectorNum"], [j])
                append!(mol.atoms.properties[num]["CycleSize"], [lastindex(subvec)])
            end
        end 
    end

    ### add further fields to atom properties
    ElemWNeighbourCount_vector = @with mol.atoms @byrow string(enumToString(:element), lastindex(neighbors(mol_graph, :number)))
    @with mol.atoms @byrow push!(:properties, 
                        "ElementWithNeighborCount" => getStringElementWithNeighborCount(:number, mol),
                        "Neighbors" => neighbors(mol_graph, :number), 
                        "AromaticityType" => mol.properties["atom_aromaticity_array"][:number],
                        "num_EWG_groups" => Int(:element) == 1 ? count_EWG(:number, mol) : -1,
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


function getStringElementWithNeighborCount(atmNum::Int64, mol::Molecule)
    return string(enumToString(mol.atoms.element[atmNum]), lastindex(neighbors(mol.properties["mol_graph"], atmNum))) 
end


function count_EWG(num::Int64, mol::Molecule)
    # only used on hydrogen atoms:
    # Electron withdrawal Atoms (EWG) according to antechamber document are: N, O, F, Cl, and Br
    # which are bound to the immediate neighbour
    # To Do: Test differences for Atom Typing, see below typical know EWG
    strong_pullers = ["Cl1", "F1", "Br1", "I1", "O1", "S1", "O2", "S2", "S3", "S4", "N2", "N3"]
    possible_indirect_pullers = ["S3", "S4", "N3"]
    acceptable_first_neighbors = ["C2", "C3", "C4", "S3", "S4", "N3", "P3", "P4", "O2", "S2"]
    elec_pullers_num = 0
    mol_graph = mol.properties["mol_graph"]
    direct_neighbor = getStringElementWithNeighborCount(neighbors(mol_graph, num)[1], mol)
    secondary_neighbors = filter(x -> !(x in neighborhood(mol_graph, num, 1)), neighborhood(mol_graph, num, 2))
    if in(direct_neighbor, acceptable_first_neighbors)
        for secNeigh in secondary_neighbors
            if getStringElementWithNeighborCount(secNeigh, mol) in strong_pullers
                elec_pullers_num += 1
            elseif getStringElementWithNeighborCount(secNeigh, mol) in possible_indirect_pullers
                tert_neighbors_elements = Vector{String}()
                for tertNeigh in filter(x -> !(x in neighborhood(mol_graph, num, 2)) && x in neighbors(mol_graph, secNeigh), 
                                                neighborhood(mol_graph, num, 3))
                    push!(tert_neighbors_elements, getStringElementWithNeighborCount(tertNeigh, mol))
                end
                if true in in(strong_pullers).(tert_neighbors_elements)
                    elec_pullers_num += 1
                end
            end
        end
    end
    return elec_pullers_num
end


function atom_conjugated_system_processor(allCycles_vec::Vector{Vector{Int64}}, mol::Molecule)
    mol_graph = mol.properties["mol_graph"]
    aromaticity_array = mol.properties["atom_aromaticity_array"]
    bond_matrix = mol.properties["weighted_graph_adj_matrix"]
    all_cycle_atoms = isempty(allCycles_vec) ? Vector{Int}() : reduce(vcat, allCycles_vec)
    conjugated_atoms_vec = Vector{Int}()
    
    filtered_bonds_df = mol.bonds[(in([Elements.C, Elements.N, Elements.O, Elements.S, Elements.P]).(mol.atoms.element[mol.bonds.a1]) .&&
                                in([Elements.C, Elements.N, Elements.O, Elements.S, Elements.P]).(mol.atoms.element[mol.bonds.a2]) .&&
                                !(true in in(aromaticity_array[mol.bonds.a1]).(["AR1","AR2"])) .&& 
                                !(true in in(aromaticity_array[mol.bonds.a2]).(["AR1","AR2"])) .&& 
                                (mol.bonds.order .== BondOrder.Double .|| mol.bonds.order .== BondOrder.Triple)), :]
    charged_atoms = filter(x -> haskey(mol.atoms[x,:properties], "Charge") && mol.atoms[x,:properties]["Charge"] != Float32(0), mol.atoms.number)
    possible_conjugated_atoms = keys(countmap(vcat(filtered_bonds_df.a1, filtered_bonds_df.a2, charged_atoms), alg=:dict))
    for atmNum in possible_conjugated_atoms
        if atmNum in charged_atoms
            push!(conjugated_atoms_vec, atmNum)
            continue
        end

        # oxygen neighbors
        oxygen_neighbors = filter(x -> mol.atoms.element[x] == Elements.O, neighbors(mol_graph, atmNum))
        oxygen_bonds = innerjoin(DataFrame(:a1 => vcat(repeat([atmNum], lastindex(oxygen_neighbors)), oxygen_neighbors), 
                                        :a2 => vcat(oxygen_neighbors, repeat([atmNum], lastindex(oxygen_neighbors)))), 
                                        filtered_bonds_df, on = [:a1, :a2])
        oxygen_bond_order_dict = countmap(oxygen_bonds.order, alg=:dict)

        if in(mol.atoms.element[atmNum], [Elements.C, Elements.N, Elements.S, Elements.P]) && nrow(oxygen_bonds) >= 2 && 
                (haskey(oxygen_bond_order_dict, BondOrder.T(2)) && oxygen_bond_order_dict[BondOrder.T(2)] >= 2)
            append!(conjugated_atoms_vec, keys(countmap(vcat(oxygen_bonds[(oxygen_bonds.order .== BondOrder.T(2)),:a1], 
                                                        oxygen_bonds[(oxygen_bonds.order .== BondOrder.T(2)),:a2]), alg=:dict)))
        elseif (in(mol.atoms.element[atmNum], [Elements.N, Elements.S, Elements.P]) &&
                (haskey(oxygen_bond_order_dict, BondOrder.T(2)) && oxygen_bond_order_dict[BondOrder.T(2)] == 1 || !haskey(oxygen_bond_order_dict, BondOrder.T(2)))) ||
                mol.atoms.element[atmNum] == Elements.C && !haskey(oxygen_bond_order_dict, BondOrder.T(2))
            direct_bonds_countmap = countmap(bond_matrix[atmNum, neighbors(mol_graph, atmNum)], alg=:dict)
            if ((haskey(direct_bonds_countmap, 2.0) && direct_bonds_countmap[2.0] >= 1) ||
                (haskey(direct_bonds_countmap, 3.0) && direct_bonds_countmap[3.0] >= 1))
                for neigh in filter(x -> !(true in in(aromaticity_array[x]).(["AR1","AR2"])), neighbors(mol_graph, atmNum))
                    neighbors_bonds_countmap = countmap(bond_matrix[neigh, 
                                                            filter(x -> !in(x, neighborhood(mol_graph, atmNum, 1)), neighbors(mol_graph, neigh))], alg=:dict)
                    if !isempty(neighbors_bonds_countmap) && 
                            in(mol.atoms.element[neigh], [Elements.C, Elements.N, Elements.O, Elements.S, Elements.P]) &&
                            ((haskey(neighbors_bonds_countmap, 2.0) && neighbors_bonds_countmap[2.0] >= 1) ||
                            (haskey(neighbors_bonds_countmap, 3.0) && neighbors_bonds_countmap[3.0] >= 1))
                        push!(conjugated_atoms_vec, atmNum)
                    end
                end
            end
        end
    end
    return conjugated_atoms_vec
end


function bondShortOrder_types(num::Int64, mol::Molecule, mol_graph::Graph, wgraph_adj_matrix::Graphs.SparseMatrix)
    BondShortVec = Vector{String}()
    bonds_dict = countmap(wgraph_adj_matrix[num, neighbors(mol_graph, num)], alg=:dict)
    for i in keys(bonds_dict)
        curr_bond_str = enumToString(BondShortOrderType(Int(i)))
        if !in("NR", mol.properties["atom_aromaticity_array"][num]) && !haskey(mol.atoms[num, :properties], "ring_non_conjugated_atom") &&
            !in("AR5", mol.properties["atom_aromaticity_array"][num]) || in(num, mol.properties["atom_conjugated_system_array"])
            push!(BondShortVec, curr_bond_str, string(bonds_dict[i], curr_bond_str))
        else
            push!(BondShortVec, uppercase(curr_bond_str), string(bonds_dict[i], uppercase(curr_bond_str)))
        end
    end
    for prop in mol.properties["atom_aromaticity_array"][num]
        push!(BondShortVec, prop)
    end
    if in(num, mol.properties["atom_conjugated_system_array"])
        push!(BondShortVec, "DL", string(countmap(mol.properties["atom_conjugated_system_array"], alg=:dict)[num], "DL"))
    elseif true in in(mol.properties["atom_aromaticity_array"][num]).(["AR1", "AR2"]) && !haskey(mol.atoms[num, :properties], "ring_non_conjugated_atom")
        push!(BondShortVec, "DL", "1DL")
    else
        push!(BondShortVec, "0DL")
    end
    return BondShortVec
end


function amide_processor(mol::Molecule, mol_graph::Graph, ElementWNeighbourCount_vector::Vector{String})
    amide_bond_vector = Vector{Tuple{Int64, Int64}}()
    nitrogen_atoms_array = mol.atoms[(mol.atoms[!, :element] .== Elements.N .&& in("NG", mol.properties["atom_aromaticity_array"][mol.atoms.number])), :number]
    for nitrogen in nitrogen_atoms_array
        for amide_neigh in neighbors(mol_graph, nitrogen)
            if lastindex(neighbors(mol_graph, amide_neigh)) > 2 && 
                in("O1", ElementWNeighbourCount_vector[neighbors(mol_graph, amide_neigh)])
                push!(amide_bond_vector, (amide_neigh,nitrogen))
            end
        end
    end
    return amide_bond_vector
end


function ClearPreprocessingMolecule!(mol::Molecule)
    mol_props_names = ["mol_graph", "adjacency_matrix", "mol_weighted_graph", 
                        "weighted_graph_adj_matrix", "chem_cycle_array", "atom_aromaticity_array", 
                        "atmprops_df", "atom_conjugated_system_array"]
    atom_props_names = ["CycleVectorNum", "CycleSize", "ElementWithNeighborCount",
                        "AromaticityType", "BondTypes", "Neighbors", "num_EWG_groups"]
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


function create_atom_preprocessing_df!(mol::Molecule)
    # Create DataFrame for better accessibility and handling of atom properties 
    
    col_names = ["AromaticityType", "ElementWithNeighborCount", "Neighbors", "num_EWG_groups", 
                 "BondTypes", "CycleSize", "CycleVectorNum"]
    num_atms = nrow(mol.atoms)
    atom_props_df = DataFrame([Vector{Vector{String}}(undef, num_atms), Vector{String}(undef, num_atms), 
                        Vector{Vector{Int64}}(undef, num_atms), Vector{Int64}(undef, num_atms), 
                        Vector{Vector{String}}(undef, num_atms), Vector{Vector{Int64}}(undef, num_atms), 
                        Vector{Vector{Int64}}(undef, num_atms)], col_names)        
    
    for atm in eachrow(mol.atoms)
        for col in intersect(col_names, keys(atm.properties))
            atom_props_df[!,col][atm.number] = atm.properties[col]
        end
    end
    return atom_props_df
end


function atom_aromaticity_type_processor(allCycles_vec::Vector{Vector{Int64}}, mol::Molecule)
    mol_graph = mol.properties["mol_graph"]
    atom_ring_class_array = Vector{Vector{String}}(undef, nrow(mol.atoms))
    reduced_vec = isempty(allCycles_vec) ? Vector{Int}() : reduce(vcat, allCycles_vec)
    for i = (1:nrow(mol.atoms))
        if i in reduced_vec
            atom_ring_class_array[i] = []
        else
            atom_ring_class_array[i] = ["NR"]
        end
    end
    for vertices_vec in allCycles_vec

        # check if is O, N, or S present in Ring vertex Vector
        ONSP_present = true in in(mol.atoms.element[vertices_vec]).([Elements.O,Elements.N,Elements.S,Elements.P])

        double_bond_to_non_ring_atom = mol.bonds[((in(vertices_vec).(mol.bonds.a1) .|| in(vertices_vec).(mol.bonds.a2)) .&&
                                                (.!in(reduced_vec).(mol.bonds.a1) .|| .!in(reduced_vec).(mol.bonds.a2)) .&&
                                                mol.bonds.order .== BondOrder.T(2)), :]
        
        if !isempty(double_bond_to_non_ring_atom)
            ring_not_conjugated_atoms = DataFrame("number" => filter(x -> in(x, vertices_vec), vcat(double_bond_to_non_ring_atom.a1, double_bond_to_non_ring_atom.a2)))
            @with innerjoin(ring_not_conjugated_atoms, mol.atoms, on = [:number]) begin
                push!.(:properties, "ring_non_conjugated_atom" => true)
            end
        end
        
        # check number of pi electrons
        subVector = copy(vertices_vec)
        shifted_subVector = vcat(subVector[2:lastindex(subVector)], subVector[1])
        ring_bonds_df = innerjoin(DataFrame(:a1 => [subVector; shifted_subVector], :a2 => [shifted_subVector; subVector]), mol.bonds, on = [:a1, :a2])
        ring_bonds_with_tripos_ar_tag_df = filter(:properties => x -> haskey(x,"TRIPOS_tag") && x["TRIPOS_tag"] == "ar", ring_bonds_df)
        ring_bonds_order_2 = haskey(countmap(ring_bonds_df.order, alg=:dict), BondOrder.T(2)) ? countmap(ring_bonds_df.order, alg=:dict)[BondOrder.T(2)] : 0
        pi_electrons = (ring_bonds_order_2 > 0 || nrow(ring_bonds_with_tripos_ar_tag_df) > 0) ? (ring_bonds_order_2 + nrow(ring_bonds_with_tripos_ar_tag_df)) * 2 : 0
        if (pi_electrons / lastindex(vertices_vec)) == 1.0
            assign_all_in_vertices_vec_aromaticity_type("AR1", vertices_vec, reduced_vec, atom_ring_class_array)
        elseif (pi_electrons / lastindex(vertices_vec)) == 0
            assign_all_in_vertices_vec_aromaticity_type("AR5", vertices_vec, reduced_vec, atom_ring_class_array)
        elseif (pi_electrons / lastindex(vertices_vec)) > 1/2 && (pi_electrons / lastindex(vertices_vec)) < 1 && 
                lastindex(vertices_vec) > 5 && isempty(double_bond_to_non_ring_atom)
            # check if Ring vertices_vec has intersections with other rings in molecule and if these are aromatic
            intersecting_atoms = filter(y -> in(y, vertices_vec), keys(countmap(filter(x -> countmap(reduced_vec)[x] > 1, reduced_vec), alg=:dict)))
            potentially_aromatic_intersecting_atoms = filter(x -> lastindex(filter(y -> y in intersecting_atoms, neighbors(mol_graph, x))) > 0, intersecting_atoms)

            check_num_aromatic_bonds = Vector{Int}()
            for potArInterAtm in potentially_aromatic_intersecting_atoms
                aromatic_intersection_bonds = mol.bonds[((mol.bonds.a1 .== potArInterAtm .|| mol.bonds.a2 .== potArInterAtm) .&&
                                                (in(filter(x -> in(neighbors(mol_graph, potArInterAtm)).(x), reduced_vec)).(mol.bonds.a1) .|| 
                                                in(filter(x -> in(neighbors(mol_graph, potArInterAtm)).(x), reduced_vec)).(mol.bonds.a2))),:order]
                ar_intersection_bonds_dict = countmap(aromatic_intersection_bonds, alg=:dict)
                push!(check_num_aromatic_bonds, haskey(ar_intersection_bonds_dict, BondOrder.T(2)) 
                                                ? ar_intersection_bonds_dict[BondOrder.T(2)]
                                                : 0)
            end
            if !isempty(check_num_aromatic_bonds) && all(in([1]).(check_num_aromatic_bonds))
                assign_all_in_vertices_vec_aromaticity_type("AR1", vertices_vec, reduced_vec, atom_ring_class_array)
            end
        elseif (pi_electrons / lastindex(vertices_vec)) > 1/2 && (pi_electrons / lastindex(vertices_vec)) <= 1 && 
                ONSP_present && lastindex(vertices_vec) > 4 && isempty(double_bond_to_non_ring_atom)
            assign_all_in_vertices_vec_aromaticity_type("AR2", vertices_vec, reduced_vec, atom_ring_class_array)
        elseif (pi_electrons / lastindex(vertices_vec)) != 0 && !isempty(double_bond_to_non_ring_atom) 
            assign_all_in_vertices_vec_aromaticity_type("AR3", vertices_vec, reduced_vec, atom_ring_class_array)
        else
            assign_all_in_vertices_vec_aromaticity_type("AR4", vertices_vec, reduced_vec, atom_ring_class_array)
        end
    end
    return atom_ring_class_array
end


function assign_all_in_vertices_vec_aromaticity_type(aromaticity::String, vertices_vec::Vector{Int}, reduced_vec::Vector{Int}, atom_ring_class_array::Vector{Vector{String}})
    for x in vertices_vec
        if !in(aromaticity, atom_ring_class_array[x])
            append!(atom_ring_class_array[x], [aromaticity], [string("RG", lastindex(vertices_vec)), string(1, "RG", lastindex(vertices_vec))])
        elseif in(aromaticity, atom_ring_class_array[x]) && 
                !isnothing(findnext(x -> x == string("RG", lastindex(vertices_vec)), atom_ring_class_array[x], 
                    !isnothing(findfirst(y -> y == aromaticity, atom_ring_class_array[x])) ? findfirst(y -> y == aromaticity, atom_ring_class_array[x]) : 1))
            count_ring_index = findnext(x -> x == string("RG", lastindex(vertices_vec)), atom_ring_class_array[x], 
                                            findfirst(y -> y == aromaticity, atom_ring_class_array[x]))+1
            curr_count = parse(Int,atom_ring_class_array[x][count_ring_index][1])
            max_count_occurences = countmap(reduced_vec, alg=:dict)[x]
            atom_ring_class_array[x][count_ring_index] = string(curr_count+1 >= max_count_occurences
                                                                ? max_count_occurences
                                                                : curr_count+1, "RG", lastindex(vertices_vec))
        elseif in(aromaticity, atom_ring_class_array[x]) && 
                isnothing(findnext(x -> x == string("RG", lastindex(vertices_vec)), atom_ring_class_array[x], 
                    !isnothing(findfirst(y -> y == aromaticity, atom_ring_class_array[x])) ? findfirst(y -> y == aromaticity, atom_ring_class_array[x]) : 1))
            insert!(atom_ring_class_array[x], findfirst(b -> b == aromaticity, atom_ring_class_array[x])+1, string("RG", lastindex(vertices_vec)))
            insert!(atom_ring_class_array[x], findfirst(b -> b == aromaticity, atom_ring_class_array[x])+2, string(1, "RG", lastindex(vertices_vec)))
        end
    end
end


function cycle_intersections(allCycles_vec::Vector{Vector{Int64}})
    inters_matrix = Matrix{Vector{Int64}}(undef, lastindex(allCycles_vec), lastindex(allCycles_vec))
    for i = (1:lastindex(allCycles_vec))
        curr1_array = copy(allCycles_vec[i])
        for j = (1:lastindex(allCycles_vec))
            curr2_array = copy(allCycles_vec[j])
            if curr1_array != curr2_array
                inters_matrix[i,j] = inters_matrix[j,i] = intersect(curr1_array, curr2_array)
            else
                inters_matrix[i,j] = inters_matrix[j,i] = []
            end
        end
    end
    return inters_matrix
end


function build_weighted_graph(mol::Molecule)
    mol_weighted_graph = SimpleWeightedGraph(nrow(mol.atoms))
    for i = (1:nrow(mol.bonds))
        add_edge!(mol_weighted_graph, mol.bonds.a1[i], mol.bonds.a2[i], Int(mol.bonds.order[i]))
    end
    return mol_weighted_graph
end


function build_graph(mol::Molecule)
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
