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
using JLD2, CodecZlib
using Logging

include("helpers/run_lDGA_dir.jl")
dir = ARGS[1]
out_fname = ARGS[2]
if !isfile(out_fname)
    jldopen(out_fname, "w"; compress = true) do f
        f["done"] = []
    end
end

jldopen(out_fname, "a+") do f_out
    done_list = f_out["done"]
    print("Processed: ", lpad(0.0, 5, " "), "%")
    i = 0
    flist = filter(x-> endswith(x,".jld2"), readdir(dir))
    for f_in_n in flist
        print("\rProcessed: ", lpad(round(100*i/length(flist),digits=1), 5, " "), "%")
        flush(stdout)
        if !(f_in_n in done_list)
            jldopen(joinpath(dir,f_in_n),"r") do f_in
                U,beta,Nk = jldopen(f_in["in_file"]) do lDGA_f
                    U = lDGA_f["mP"].U
                    beta = lDGA_f["mP"].β
                    Nk = lDGA_f["Nk"]
                    if !haskey(f_out, "$U@$beta")
                        f_out["$U@$beta/beta"] = beta
                        f_out["$U@$beta/U"] = U
                    end
                    f_out["$U@$beta/$Nk/Nk"] = Nk
                    f_out["$U@$beta/$Nk/χ_sp"] = lDGA_f["nlQ_sp"].χ
                    f_out["$U@$beta/$Nk/χ_ch"] = lDGA_f["nlQ_ch"].χ
                    f_out["$U@$beta/$Nk/lDGAlog"] = lDGA_f["log"]
                    f_out["$U@$beta/$Nk/χ_sp_usable"] = lDGA_f["nlQ_sp"].usable_ω
                    f_out["$U@$beta/$Nk/χ_ch_usable"] = lDGA_f["nlQ_ch"].usable_ω
                    f_out["$U@$beta/$Nk/λsp_old"] = lDGA_f["λsp_old"]
                    U,beta,Nk
                end
                f_out["$U@$beta/$Nk/res"] = f_in["res"]
                f_out["$U@$beta/$Nk/λsp"] = f_in["λsp_new"]
                f_out["$U@$beta/$Nk/λch"] = f_in["λch_new"]
                f_out["$U@$beta/$Nk/Epot_DMFT_0"] = f_in["Epot_DMFT_0"]
                f_out["$U@$beta/$Nk/Epot_DMFT_1"] = f_in["Epot_DMFT_1"]
                f_out["$U@$beta/$Nk/Epot_DMFT_2"] = f_in["Epot_DMFT_2"]
                f_out["$U@$beta/$Nk/Epot_lDGA_λsp_1"] = f_in["Epot_lDGA_λsp_1"]
                f_out["$U@$beta/$Nk/Epot_lDGA_λsp_2"] = f_in["Epot_lDGA_λsp_2"]
                f_out["$U@$beta/$Nk/Epot_lDGA_λspch_1"] = f_in["Epot_lDGA_λspch_1"]
                f_out["$U@$beta/$Nk/Epot_lDGA_λspch_2"] = f_in["Epot_lDGA_λspch_2"]
                f_out["$U@$beta/$Nk/Ekin_DMFT_0"] = f_in["Ekin_DMFT_0"]
                f_out["$U@$beta/$Nk/Ekin_DMFT_1"] = f_in["Ekin_DMFT_1"]
                f_out["$U@$beta/$Nk/Ekin_lDGA_λsp_1"] = f_in["Ekin_lDGA_λsp_1"]
                f_out["$U@$beta/$Nk/Ekin_lDGA_λspch_1"] = f_in["Ekin_lDGA_λspch_1"]
                f_out["$U@$beta/$Nk/log"] = f_in["log"]
                push!(done_list,f_in_n)
                Base.delete!(f_out, "done")
                f_out["done"] = done_list
            end
        end
        i += 1
    end
    println("\rDone.")
end
