"""
    Description:
Reads various results from ladderDGA.jl results from files starting with ARG2 in subdirectory of ARG1. Scans for _kNUMBER_ in file name and reads data for all of these seperately. Finally writes combined jld2 results for all beta,U and Nk to ARG3.jld2

    Example:
julia PD_gather_ext_l.jl /scratch/projects/hhp00048/lDGA/PD/data lDGA_ res
"""

using Pkg
Pkg.activate("/scratch/projects/hhp00048/codes/LadderDGA.jl",io=devnull)
using JLD2
using HDF5
using TOML
using LadderDGA
using DataFrames
using NLsolve
include("helpers/new_lambda_helpers.jl")


dir = ARGS[1]
fname_pre = ARGS[2]
out_fname = ARGS[3]

df = DataFrame(β = Float64[], U = Float64[], Nk = Int[], tc=Symbol[],
               χsp = Matrix{ComplexF64}[], χch = Matrix{ComplexF64}[], 
               λch_range = Vector{Float64}[], spOfch = Vector{Float64}[], λsp_of_λch_res = Matrix{Float64}[], 
               λsp = Float64[], λnew_sp = Float64[], λnew_ch = Float64[], λnew_converged = Bool[],
               λsp_EPot_intern = Float64[], λnew_sp_EPot_intern = Float64[], λnew_ch_EPot_intern = Float64[],
               EPot_direct_DMFT = Float64[], EPot_GS_DMFT = Float64[],
               EPot_chi_λ0 = Float64[], EPot_GS_λ0 = Float64[],
               EPot_chi_λsp = Float64[], EPot_GS_λsp = Float64[],
               EPot_chi_λspch = Float64[], EPot_GS_λspch = Float64[],
               config = String[], log = String[])

dpw = displaysize(stdout)
println("Walking through $dir")
for (root, dirs, files) in walkdir(dir)
    if "config.toml" in files
        flist = filter(x->startswith(x,fname_pre), files)
        print("\r$(repeat(" ", dpw[2]))")
        for file in flist
            m = match(r"_k(?<kn>\d+)_",file)
            if m !== nothing
                kn = parse(Int,m[:kn])     
                #TODO: DO NOT HARDCODE LATTICE TYPE!!!
                jldopen(joinpath(root, file), "r") do f
                    rr = E_pot_test(f["bubble"], f["impQ_sp"], f["impQ_ch"], f["nlQ_sp"], f["nlQ_ch"], f["Sigma_loc"], kn, f["gLoc"], f["gLoc_fft"], f["Sigma_DMFT"], f["FUpDo"], f["mP"], f["sP"])
                    res_nls = f["λnew_nls"]
                    r = [f["mP"].β, f["mP"].U, f["kG"].Ns, f["sP"].tc_type_f,
                         f["nlQ_sp"].χ, f["nlQ_ch"].χ,
                         collect(f["λch_range"]), f["spOfch"], f["λsp_of_λch_res"], f["λsp_old"], 
                         res_nls.zero..., res_nls.f_converged, 
                         rr[4:(end-2)]...,
                         f["config"], f["log"]]
                    push!(df, r)
                end
            end
        end
    end
end

jldopen(out_fname,"a+") do f_out
    f_out["df"] = df
end
