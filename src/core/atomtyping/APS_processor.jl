# APS is the atomic_property string from the DEF file format.
# In dataframes of imported DEF files there is a column "atomic_property"

export APS_processor

function APS_processor(colstring::String, atmprops_df::DataFrameRow)
    colstring = colstring[2:lastindex(colstring)-1]
    and_count = count(==(','), colstring)
    or_count = count(==('.'), colstring)
    and_list = split(colstring, ',')
    and_expr_conditionals = Expr(:&&)

    if !all(in(atmprops_df.BondTypes).(and_list))
        for andItem in and_list
            or_list = split(andItem, '.')
            or_expr_conditionals = Expr(:||)
            for orItem in or_list
                push!(or_expr_conditionals.args, :(in($(atmprops_df.BondTypes)).($orItem)))
            end
            if eval(or_expr_conditionals) == false
                return false
            else
                push!(and_expr_conditionals.args, eval(or_expr_conditionals)) 
            end
        end
    else
        return true
    end
    return eval(and_expr_conditionals)
end