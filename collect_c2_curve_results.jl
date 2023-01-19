"""
    Description:
Read through directory tree given by ARG1, containing results from `run_c2_curves.jl`.
Combines results to file, given by ARG2.

    Example:
julia collect_c2_curve_results.jl res res.jld2
"""

using Pkg
Pkg.activate("/scratch/projects/hhp00048/codes/LadderDGA.jl")
using LadderDGA
using Dispersions
using TOML
using JLD2
using Logging
using DataFrames

include("helpers/run_lDGA_dir.jl")
dir = ARGS[1]
fname_pre = ARGS[2]
out_fname = ARGS[2]
if !isfile(out_fname)
    jldopen(out_fname, "w"; compress = true) do f
        f["done"] = []
    end
end

df = DataFrame(β = Float64[], U = Float64[], ndens = Float64[],
               λsp_old = Float64[], λsp = Float64[], λch = Float64[],
               λsp_sc = Float64[], λch_sc = Float64[],
               χAF = Float64[], χAF_m = Float64[], χAF_dm = Float64[],
               χAF_dm_sc = Float64[])


println("Walking through $dir")
for (root, dirs, files) in walkdir(dir)
    flist = filter(x->startswith(x,fname_pre), files)
    flist = filter(x->endswith(x,".jld2"), flist)
    for file in flist
        println("file $(joinpath(root, file))")
        flush(stdout)
        jldopen(joinpath(root, file), "r") do f
            U = f["mP"].U
            beta = f["mP"].β
            sP = f["sP"]
            mP = f["mP"]
            χAF = real(f["χsp"].data[end,151])
            χAF_m = 1/ (1 / χAF + f["λsp_old"])
            χAF_dm = 1/ (1 / χAF + f["λspch"][1])
            χAF_dm_sc = 1/ (1 / χAF + f["λspch_sc"][1])

            row = [beta, U, mP.n, f["λsp_old"], f["λspch"][1], f["λspch"][2], f["λspch_sc"][1], f["λspch_sc"][2], χAF, χAF_m, χAF_dm, χAF_dm_sc]
            push!(df, row)
        end
    end
end

println(df)
