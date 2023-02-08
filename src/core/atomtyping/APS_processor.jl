# APS is the atomic_property string from the DEF file format.
# In dataframes of imported DEF files there is a column "atomic_property"

export APS_processor

function APS_processor(colstring::String, atmprops_df::DataFrameRow)
    if colstring[1] == '[' && colstring[lastindex(colstring)] == ']'
        colstring = colstring[2:lastindex(colstring)-1]
    end
    
    and_list = split(colstring, ',')
    and_Bool = true

    if !all(in(atmprops_df.BondTypes).(and_list))
        for andItem in and_list
            or_list = split(andItem, '.')
            or_Bool = false
            for orItem in or_list
                or_Bool = in(atmprops_df.BondTypes).(orItem)
                if or_Bool == true
                    break
                end
            end
            if or_Bool == false
                return false
            end
        end
    else
        return true
    end
    return and_Bool
end