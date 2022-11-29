# APS is the atomic_property string from the DEF file format.
# In dataframes of imported DEF files there is a column "atomic_property"


function APS_processor(colstring::String, atmprops_df::DataFrameRow, splitchar::Char)
    and_count = count(==(','), colstring)
    or_count = count(==('.'), colstring)
    and_list = split(colstring, ',')
    and_expr_conditionals = :(:filler)

    if !all(in(atmprops_df.BondTypes).(and_list))
        for num in 1:and_count
            if num > 1
                push!(and_expr_conditionals.args, :(&&), :filler)
            end
        end
        for (num, item) in enumerate(and_list)
            or_list = split(item, '.')
            or_expr_conditionals = :(:filler)
            for (j, orItem) in enumerate(or_list)
                if j == 1
                   or_expr_conditionals.args[1] = :(in(atmprops_df.BondTypes).(orItem)) 
                else
                    push!(or_expr_conditionals.args, :(||), :(in(atmprops_df.BondTypes).(orItem)))
                end
            end
            and_expr_conditionals.args[2*(num-1)+1] = eval(or_expr_conditionals) 
        end
    else
        return true
    end
    return eval(and_expr_conditionals)
end