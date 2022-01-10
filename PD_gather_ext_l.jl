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


dir = ARGS[1]
fname_pre = ARGS[2]
out_fname = ARGS[3]


jldopen(out_fname*".jld2", "a+") do f_out
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
                    jldopen(joinpath(root, file)) do f
                        U = f["mP"].U
                        beta = f["mP"].β
                        if !haskey(f_out, "$U@$beta")
                            f_out["$U@$beta/beta"] = beta
                            f_out["$U@$beta/U"] = U
                        end
                        println("file $(joinpath(root, file)) with U=$U, β=$beta and λsp=$(f["λsp_old"])")

                        f_out["$U@$beta/$kn/Nk"] = kn
                        f_out["$U@$beta/$kn/chi_sp_ppp"] = f["nlQ_sp"].χ[end,:]
                        f_out["$U@$beta/$kn/chi_ch_ppp"] = f["nlQ_ch"].χ[end,:]
                        f_out["$U@$beta/$kn/chi_ch_000"] = f["nlQ_ch"].χ[1,:]
                        f_out["$U@$beta/$kn/lambda_sp_old"] = f["λsp_old"]
                        f_out["$U@$beta/$kn/λch_range"] = f["λch_range"]
                        f_out["$U@$beta/$kn/spOfch"] = f["spOfch"]
                        f_out["$U@$beta/$kn/λsp_of_λch_res"] = f["λsp_of_λch_res"]
                    end
                end
            end
        end
    end
end
