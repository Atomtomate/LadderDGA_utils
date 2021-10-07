function filter_usable_λsp_of_λch(λch_range, inp; max_λsp=Inf)
    tmp = deepcopy(inp)
    tmp[isnan.(tmp)] .= 0.0
    tmp[tmp .> max_λsp] .= 0.0
    findmax(tmp)[2]:length(λch_range)
end
