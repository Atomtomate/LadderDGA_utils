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

include("../helpers/run_lDGA_dir.jl")
dir = ARGS[1]
fname_pre = ARGS[2]
out_fname = ARGS[3]


jldopen(out_fname*".jld2", "a+") do f_out
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
                ndens = mP.n
                key = "$U@$beta@$ndens"
                f_out[key*"/U"] = U
                f_out[key*"/β"] = beta
                f_out[key*"/ndens"] = ndens
                f_out[key*"/sP"] = sP
                f_out[key*"/mP"] = mP
                f_out[key*"/Σ_loc"] = f["Σ_loc"]
                f_out[key*"/Σ_ladder_m"] = f["Σ_ladder_m"]
                f_out[key*"/Σ_ladder_dm"] = f["Σ_ladder_dm"]
                f_out[key*"/Σ_ladder_m_sc"] = f["Σ_ladder_m_sc"]
                f_out[key*"/Σ_ladder_dm_sc"] = f["Σ_ladder_dm_sc"]
                f_out[key*"/Σ_ladder_m_tsc"] = f["Σ_ladder_m_tsc"]
                f_out[key*"/Σ_ladder_dm_tsc"] = f["Σ_ladder_dm_tsc"]
                f_out[key*"/converged_m"] = f["converged_m"]
                f_out[key*"/converged_m_sc"] = f["converged_m_sc"]
                f_out[key*"/converged_dm"] = f["converged_dm"]
                f_out[key*"/converged_dm_sc"] = f["converged_dm_sc"]
                f_out[key*"/converged_m_tsc"] = f["converged_m_tsc"]
                f_out[key*"/converged_dm_tsc"] = f["converged_dm_tsc"]
                f_out[key*"/ef_m"] = f["ef_m"]
                f_out[key*"/ef_dm"] = f["ef_dm"]
                f_out[key*"/ef_dm_sc"] = f["ef_dm_sc"]
                f_out[key*"/ef_dm_tsc"] = f["ef_dm_tsc"]
            end
        end
    end
end
