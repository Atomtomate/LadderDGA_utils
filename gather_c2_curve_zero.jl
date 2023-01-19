"""
    Description:
Reads various results from ladderDGA.jl results from files starting with ARG2 in subdirectory of ARG1. Scans for _kNUMBER_ in file name and reads data for all of these seperately. Finally writes combined jld2 results for all beta,U and Nk to ARG3.jld2

    Example:
julia PD_gather_kConv.jl /scratch/projects/hhp00048/lDGA/PD/data lDGA_ res
"""

using Pkg
Pkg.activate("/scratch/projects/hhp00048/codes/LadderDGA.jl",io=devnull)
using JLD2
using TOML
using LadderDGA


dir = ARGS[1]
fname_pre = ARGS[2]
out_fname = ARGS[3]
error_list = []


jldopen(out_fname*".jld2", "a+") do f_out
    dpw = displaysize(stdout)
    println("Walking through $dir")
    for (root, dirs, files) in walkdir(dir)
        flist = filter(x->startswith(x,fname_pre), files)
        flist = filter(x->endswith(x,".jld2"), flist)
        #print("\r$(repeat(" ", dpw[2]))")
        for file in flist
            println("file $(joinpath(root, file))")
            flush(stdout)
            jldopen(joinpath(root, file), "r") do f
                cfg = TOML.parse(f["config"]) 
                U = cfg["Model"]["U"]
                β = cfg["Model"]["beta"]
                ndens = cfg["Model"]["nden"]
                root_key = "$U@$β@$ndens"
                if haskey(f_out, root_key)
                    println("WARNING: duplicate key $root_key for file $root/$file")
                    return
                end
                f_out["λsp_old"] = f["λsp_old"] 
                f_out["λsp_new"] = f["λsp_new"]
                f_out["λch_new"] = f["λch_new"]
                f_out["check_new"] = f["check_new"]
                f_out["λsp_new_sc"] = f["λsp_new_sc"]
                f_out["λch_new_sc"] = f["λch_new_sc"]
            end
        end
    end
end

println("ERROR LIST: \n", error_list)
