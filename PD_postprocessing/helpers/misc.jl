function dbg_fix_bad_Nk_group(f, target::Int)
    for k in keys(f)
        for ki in keys(f[k])
            println("in $k/$ki with tryparse $(tryparse(Int, ki)) and target $(tryparse(Int, ki) != nothing && parse(Int, ki) != target)")
            if ((tryparse(Int, ki) != nothing) && (parse(Int, ki) != target))
                println("trying to move $(k*"/"*ki) to $(k)/$target" )
                f["$k/$(target)"] = deepcopy(f[k*"/"*ki])
                delete!(f, k*"/"*ki)
            end
        end
    end
end
