
using DataFrames

export load_atomtyping_DEF

function load_atomtyping_DEF(mapfile::AbstractString)
    ATD_lines_vec = filter(x -> lastindex(x) > 3 && x[1:3] == "ATD", readlines(mapfile))
    df = DataFrame([Vector{Union{String, Nothing}}(undef, lastindex(ATD_lines_vec)), Vector{Union{String, Nothing}}(undef, lastindex(ATD_lines_vec)), 
                    Vector{Union{Int, Nothing}}(undef, lastindex(ATD_lines_vec)), Vector{Union{Int, Nothing}}(undef, lastindex(ATD_lines_vec)), 
                    Vector{Union{Int, Nothing}}(undef, lastindex(ATD_lines_vec)), Vector{Union{Int, Nothing}}(undef, lastindex(ATD_lines_vec)), 
                    Vector{Union{String, Nothing}}(undef, lastindex(ATD_lines_vec)), Vector{Union{String, Nothing}}(undef, lastindex(ATD_lines_vec))],
                    ["type_name", "residue_names", "atomic_number", "num_neighbors",
                    "num_H_bonds", "electron_withdrawal_groups", "atomic_property", "CES"])

    for (line_num, line) in enumerate(ATD_lines_vec)
        line_vec = line_in_def_df(line)
        df.type_name[line_num] = line_vec[1]
        df.residue_names[line_num] = line_vec[2]
        df.atomic_number[line_num] = line_vec[3]
        df.num_neighbors[line_num] = line_vec[4]
        df.num_H_bonds[line_num] = line_vec[5]
        df.electron_withdrawal_groups[line_num] = line_vec[6]
        df.atomic_property[line_num] = line_vec[7]
        df.CES[line_num] = line_vec[8]
    end
    return df
end

function line_in_def_df(line::AbstractString)
    row_data_list = Vector{Union{AbstractString, Int64, Nothing}}()
    line_data = split(line)
    for i = (2:9)
        if (!isassigned(line_data, i) || in(line_data[i], ["&", "*"]))
            append!(row_data_list, [nothing])
        elseif isnumeric(line_data[i][1])
            append!(row_data_list, [parse(Int, line_data[i])]) 
        else
            append!(row_data_list, [line_data[i]])
        end
    end
    return row_data_list
end


