
using DataFrames

export load_atomtyping_DEF

function load_atomtyping_DEF(mapfile::AbstractString)
    df = DataFrame([Vector{Union{String, Nothing}}(), Vector{Union{String, Nothing}}(), 
                    Vector{Union{Int, Nothing}}(), Vector{Union{Int, Nothing}}(), Vector{Union{Int, Nothing}}(), Vector{Union{Int, Nothing}}(), 
                    Vector{Union{String, Nothing}}(), Vector{Union{String, Nothing}}()],
                    ["type_name", "residue_names", "atomic_number", "num_neighbors",
                    "num_H_bonds", "electron_withdrawal_groups", "atomic_property", "CES"])
    
    for line in filter(x -> lastindex(x) > 3 && x[1:3] == "ATD", readlines(mapfile))
        push!(df, lines_for_df(line))
    end
    
    return df
end

function lines_for_df(line::AbstractString)
    row_data_list = Vector{Union{AbstractString, Int64, Nothing}}()
    line_data = split(line)
    for i = (2:9)
        if (!isassigned(line_data, i) || in(["&", "*"]).(line_data[i]))
            append!(row_data_list, [nothing])
        elseif isnumeric(line_data[i][1])
            append!(row_data_list, [parse(Int, line_data[i])]) 
        else
            append!(row_data_list, [line_data[i]])
        end
    end
    return row_data_list
end


