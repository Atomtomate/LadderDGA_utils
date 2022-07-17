"""
    Description:
Reads various results from ladderDGA.jl results from files starting with ARG2 in subdirectory of ARG1. Scans for _kNUMBER_ in file name and reads data for all of these seperately. Finally writes combined jld2 results for all beta,U and Nk to ARG3.jld2

    Example:
julia PD_gather.jl /scratch/projects/hhp00048/lDGA/PD/data lDGA_ res
"""

using Pkg
Pkg.activate("/scratch/projects/hhp00048/codes/LadderDGA.jl",io=devnull)
using JLD2
using TOML
using LadderDGA


dir = ARGS[1]
fname_pre = ARGS[2]
out_fname = ARGS[3]


jldopen(out_fname*".jld2", "a+") do f_out
    dpw = displaysize(stdout)
    println("Walking through $dir")
    for (root, dirs, files) in walkdir(dir)
        if "config.toml" in files
            flist = filter(x->startswith(x,fname_pre), files)
            flist = filter(x->endswith(x,".jld2"), flist)
            #print("\r$(repeat(" ", dpw[2]))")
            for file in flist
                println("file $(joinpath(root, file))")
                flush(stdout)
                jldopen(joinpath(root, file), "r") do f
                    if !haskey(f, "kG")
                        println("WARNING: Skipping corrupt file: $file")
                        return
                    end
                    kn = f["kG"].Ns
                    U = f["mP"].U
                    beta = f["mP"].β
                    root_key = "$U@$beta"
                    if !haskey(f_out, root_key)
                        f_out[root_key*"/beta"] = beta
                        f_out[root_key*"/U"] = U
                    else
                        if haskey(f_out, root_key*"/$kn")
                            println("WARNING: Duplicate key $root_key/$(kn)! Skipping gather")
                            return
                        end
                    end
                    nh = ceil(Int,size(f["nlQ_sp"].χ, 2)/2)
                    println("   ---> U=$U, β=$beta and λsp=$(f["λsp"])")

                    f_out[root_key*"/$kn/description"] = f["Description"]
                    f_out[root_key*"/$kn/config"] = f["config"]
                    f_out[root_key*"/$kn/Nk"] = kn
                    f_out[root_key*"/$kn/Σ_ladder_DMFT"] = f["Σ_ladder_DMFT"]
                    f_out[root_key*"/$kn/Σ_ladder_λsp"] = f["Σ_ladder_λsp"]
                    f_out[root_key*"/$kn/Σ_ladder_λspch"] = f["Σ_ladder_λspch"]
                    f_out[root_key*"/$kn/Σ_DMFT"] = f["Sigma_DMFT"]
                    f_out[root_key*"/$kn/λsp"] = f["λsp"]
                    f_out[root_key*"/$kn/λspch"] = f["λspch"]

                    f_out[root_key*"/$kn/χAF_DMFT"] = haskey(f, "χAF_DMFT") ? f["χAF_DMFT"] : real(1 / (1 / f["nlQ_sp"].χ[end,nh]))
                    f_out[root_key*"/$kn/χAF_λsp"] = haskey(f, "χAF_λsp") ? f["χAF_λsp"] : real(1 / (1 / f["nlQ_sp"].χ[end,nh] + f["λsp"]))

                    f_out[root_key*"/$kn/χAF_λspch"] = if haskey(f, "χAF_λspch")
                        f["χAF_λspch"]
                    else
                    if typeof(f["λspch"]) === Array
                        real(1 / (1 / f["nlQ_sp"].χ[end,nh] + f["λspch"][1]))
                    else
                        0.0
                    end
                    end
                    f_out[root_key*"/$kn/log"] = f["log"]
                    f_out[root_key*"/$kn/c2_curve"] = f["c2_curve"]
                    f_out[root_key*"/$kn/E_kin_ED"] = f["E_kin_ED"]
                    f_out[root_key*"/$kn/E_pot_ED"] = f["E_pot_ED"]
                    f_out[root_key*"/$kn/E_kin_DMFT_1"] = f["E_kin_DMFT_1"]
                    f_out[root_key*"/$kn/E_pot_DMFT_1"] = f["E_pot_DMFT_1"]
                    f_out[root_key*"/$kn/E_pot_DMFT_2"] = f["E_pot_DMFT_2"]
                    f_out[root_key*"/$kn/E_kin_λsp_1"] = f["E_kin_λsp_1"]
                    f_out[root_key*"/$kn/E_pot_λsp_1"] = f["E_pot_λsp_1"]
                    f_out[root_key*"/$kn/E_pot_λsp_2"] = f["E_pot_λsp_2"]
                    f_out[root_key*"/$kn/E_kin_λspch_1"] = f["E_kin_λspch_1"]
                    f_out[root_key*"/$kn/E_pot_λspch_1"] = f["E_pot_λspch_1"]
                    f_out[root_key*"/$kn/E_pot_λspch_2"] = f["E_pot_λspch_2"]
                    f_out[root_key*"/$kn/c2_λspch"] = f["c2_λspch"]
                end
            end
        end
    end
end
