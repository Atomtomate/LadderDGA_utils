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
                U = f["mP"].U
                β = f["mP"].β
                root_key = "$U@$β"
                if haskey(f_out, root_key)
                    println("WARNING: duplicate key $root_key for file $root/$file")
                    return
                end
                f_out[root_key*"/E_kin_ED"] = f["E_kin_ED"]
                f_out[root_key*"/E_pot_ED"] = f["E_pot_ED"]
                f_out[root_key*"/config"] = f["config"]
                f_out[root_key*"/β"] = β
                f_out[root_key*"/U"] = U
                kn_list = filter(x-> x != nothing, tryparse.(Int, keys(f)))
                if isempty(kn_list)
                    println("WARNING: U/β = $U/$β, no data found!")
                    push!(error_list, (U,β,:no_data))
                    return
                end
                kn = maximum(kn_list)

                f_out[root_key*"/Nk"] = kn
                f_out[root_key*"/λsp"] = f["$kn/λsp"]
                f_out[root_key*"/λspch"] = f["$kn/λspch"]

                f_out[root_key*"/χAF_DMFT"] = f["$kn/χAF_DMFT"]
                f_out[root_key*"/χAF_λsp"] = f["$kn/χAF_λsp"]
                f_out[root_key*"/χAF_λspch"] = f["$kn/χAF_λspch"]
                f_out[root_key*"/log"] = f["$kn/log"]
                f_out[root_key*"/E_pot_DMFT_2"] = f["$kn/E_pot_DMFT_2"]
                f_out[root_key*"/E_kin_λsp_1"] = f["$kn/E_kin_λsp_1"]
                f_out[root_key*"/E_pot_λsp_1"] = f["$kn/E_pot_λsp_1"]
                f_out[root_key*"/E_pot_λsp_2"] = f["$kn/E_pot_λsp_2"]
                f_out[root_key*"/E_kin_λspch_1"] = f["$kn/E_kin_λspch_1"]
                f_out[root_key*"/E_pot_λspch_1"] = f["$kn/E_pot_λspch_1"]
                f_out[root_key*"/E_pot_λspch_2"] = f["$kn/E_pot_λspch_2"]
                f_out[root_key*"/conv_error"] = f["$kn/conv_error"]
                f_out[root_key*"/conv"] = f["$kn/conv"]
                if f["$kn/conv_error"] || !f["$kn/conv"]
                    println("WARNING: U/β = $U/$β did not converge!")
                    push!(error_list, (U,β,:not_converged))
                end
                if kn <= 10
                    push!(error_list, (U,β,:only_k10))
                    println("WARNING: U/β = $U/$β, max KN is $kn")
                end
            end
        end
    end
end

println("ERROR LIST: \n", error_list)
