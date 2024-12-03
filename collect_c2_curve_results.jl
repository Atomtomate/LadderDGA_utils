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

include(abspath(joinpath(@__DIR__,"../helpers/run_lDGA_dir.jl")))
dir = ARGS[1]
fname_pre = ARGS[2]
out_fname = ARGS[2]
if !isfile(out_fname)
    jldopen(out_fname, "w"; compress = true) do f
        f["done"] = []
    end
end



df = DataFrame(β = Float64[], U = Float64[], ndens = Float64[],
               μ_DMFT = Float64[], μ_m = Float64[], μ_m_sc = Float64[], μ_dm = Float64[], μ_dm_sc = Float64[], μ_m_tsc = Float64[], μ_dm_tsc = Float64[],
               λ_m = Float64[], λ_m_sc = Float64[], λ_dm_m = Float64[], λ_dm_d = Float64[], λ_dm_sc_m = Float64[], λ_dm_sc_d = Float64[], λ_m_tsc = Float64[], λ_dm_tsc_m = Float64[], λ_dm_tsc_d = Float64[],
               E_kin_DMFT = Float64[], E_kin_m = Float64[], E_kin_m_sc = Float64[], E_kin_dm = Float64[], E_kin_dm_sc = Float64[], E_kin_m_tsc = Float64[], E_kin_dm_tsc = Float64[],
               E_pot_DMFT = Float64[], E_pot_m = Float64[], E_pot_m_sc = Float64[], E_pot_dm = Float64[], E_pot_dm_sc = Float64[], E_pot_m_tsc = Float64[], E_pot_dm_tsc = Float64[],
               χ0inv_0_DMFT = Float64[], χ0inv_0_m = Float64[], χ0inv_0_m_sc = Float64[], χ0inv_0_dm = Float64[], χ0inv_0_dm_sc = Float64[], χ0inv_0_m_tsc = Float64[], χ0inv_0_dm_tsc = Float64[],
               χ0inv_π_DMFT = Float64[], χ0inv_π_m = Float64[], χ0inv_π_m_sc = Float64[], χ0inv_π_dm = Float64[], χ0inv_π_dm_sc = Float64[], χ0inv_π_m_tsc = Float64[], χ0inv_π_dm_tsc = Float64[],
               converged_m = Bool[], converged_m_sc = Bool[], converged_dm = Bool[], converged_dm_sc = Bool[], converged_m_tsc = Bool[], converged_dm_tsc = Bool[],
               χAF_DMFT = Float64[], χAF_m = Float64[], χAF_dm = Float64[], χAF_dm_sc = Float64[], χAF_m_tsc = Float64[], χAF_dm_tsc = Float64[])


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
            nh = ceil(Int,size(f["χsp"].data,2)/2) 
            χAF = real(f["χsp"].data[end,nh])
            χAF_m = 1/ (1 / χAF + f["λm"])
            χAF_dm = 1/ (1 / χAF + f["λdm"][1])
            χAF_dm_sc = 1/ (1 / χAF + f["λdm_sc"][1])
            χAF_m_tsc = 1/ (1 / χAF + f["λm_tsc"][1])
            χAF_dm_tsc = 1/ (1 / χAF + f["λdm_tsc"][1])

            row = [beta, U, mP.n, 
                   f["μ_dmft"], f["μ_m"], f["μ_m_sc"], f["μ_dm"], f["μ_dm_sc"], f["μ_m_tsc"], f["μ_dm_tsc"],
                   f["λm"], f["λdm"][1], f["λdm"][2], f["λm_sc"], f["λdm_sc"][1], f["λdm_sc"][2], f["λm_tsc"], f["λdm_tsc"][1], f["λdm_tsc"][2],
                   f["E_kin_DMFT"], f["E_kin_m"], f["E_kin_m_sc"], f["E_kin_dm"], f["E_kin_dm_sc"], f["E_kin_m_tsc"], f["E_kin_dm_tsc"],
                   f["E_pot_DMFT"], f["E_pot_m"], f["E_pot_m_sc"], f["E_pot_dm"], f["E_pot_dm_sc"], f["E_pot_m_tsc"], f["E_pot_dm_tsc"],
                   real.([f["χ0_inv_DMFT_0"], f["χ0_inv_m_0"], f["χ0_inv_m_sc_0"], f["χ0_inv_dm_0"], f["χ0_inv_dm_sc_0"], f["χ0_inv_m_tsc_0"], f["χ0_inv_dm_tsc_0"],
                          f["χ0_inv_DMFT_π"], f["χ0_inv_m_π"], f["χ0_inv_m_sc_π"], f["χ0_inv_dm_π"], f["χ0_inv_dm_sc_π"], f["χ0_inv_m_tsc_π"], f["χ0_inv_dm_tsc_π"]])...,
                   f["converged_m"], f["converged_m_sc"], f["converged_dm"], f["converged_dm_sc"], f["converged_m_tsc"], f["converged_dm_tsc"],
                   χAF, χAF_m, χAF_dm, χAF_dm_sc, χAF_m_tsc, χAF_dm_tsc
                  ]
            push!(df, row)
        end
    end
end

println(df)
