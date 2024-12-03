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
using CSV

include("../helpers/run_lDGA_dir.jl")
dir = ARGS[1]
fname_pre = ARGS[2]
out_fname = ARGS[3]

df = DataFrame(β = Float64[], U = Float64[], ndens = Float64[], Nk = Int[], μ_dmft = Float64[], μ_sc = Float64[],
               λsp_old = Float64[], λsp = Float64[], λch = Float64[],
               λsp_sc = Float64[], λch_sc = Float64[],
               χAF = Float64[], χAF_m = Float64[], χAF_dm = Float64[],χAF_dm_sc = Float64[],
               χd0 = Float64[], χd0_m = Float64[], χd0_dm = Float64[],χd0_dm_sc = Float64[],
               χdπ = Float64[], χdπ_m = Float64[], χdπ_dm = Float64[],χdπ_dm_sc = Float64[],
               χ0_i_0_dmft = Float64[], χ0_i_π_dmft = Float64[], χ0_i_0 = Float64[], χ0_i_π = Float64[] 
               )


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
            nh = ceil(Int, size(f["χsp"].data,2)/2)
            χAF = real(f["χsp"].data[end,nh])
            χAF_m = 1/ (1 / χAF + f["λsp_old"])
            χAF_dm = 1/ (1 / χAF + f["λspch"][1])
            χAF_dm_sc = 1/ (1 / χAF + f["λspch_sc"][1])
            χd0 = real(f["χch"].data[1,nh])
            χd0_m = 1/ (1 / χd0)
            χd0_dm = 1/ (1 / χd0 + f["λspch"][2])
            χd0_dm_sc = 1/ (1 / χd0 + f["λspch_sc"][2])
            χdπ = real(f["χch"].data[end,nh])
            χdπ_m = 1/ (1 / χdπ)
            χdπ_dm = 1/ (1 / χdπ + f["λspch"][2])
            χdπ_dm_sc = 1/ (1 / χdπ + f["λspch_sc"][2])

            row = [beta, U, mP.n, size(f["χsp"].data,1), mP.μ, f["μsc"], f["λsp_old"], f["λspch"][1], f["λspch"][2], f["λspch_sc"][1], f["λspch_sc"][2], χAF, χAF_m, χAF_dm, χAF_dm_sc, χd0, χd0_m, χd0_dm, χd0_dm_sc, χdπ, χdπ_m, χdπ_dm, χdπ_dm_sc, real(f["χ0_i_0_DMFT"]), real(f["χ0_i_π_DMFT"]), real(f["χ0_i_0"]), real(f["χ0_i_π"])]
            push!(df, row)
        end
    end
end

println(df)
CSV.write(out_fname, df)
