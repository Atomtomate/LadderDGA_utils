"""
    Description:
Read through directory tree given by ARG1 and search for chi_lDGA output starting with ARG2. 
Calculate condition for lambda_ch as function of lambda_ch (root of this function is the lDGA 
value of lambda_ch).
Results are written to ARG3.
ARG4 is the kIteration in the config files.
ARG5 is the string identifying the tail correction type.

    Example:
julia run_c2_curves.jl /scratch/projects/hhp00048/lDGA/PD/data lDGA_ res
"""

using Distributed
@everywhere using Pkg
@everywhere Pkg.activate("/scratch/projects/hhp00048/codes/LadderDGA.jl",io=devnull)
@everywhere using LadderDGA
@everywhere using JLD2
using DataFrames

include("helpers/run_lDGA_dir.jl")
@everywhere include("helpers/new_lambda_helpers.jl")
@everywhere include("new_lambda_analysis.jl")

@everywhere dir = $(ARGS[1])
@everywhere fname_pre = $(ARGS[2])
@everywhere out_fname = $(ARGS[3])
@everywhere kIt = $(ARGS[4])
@everywhere tcType = $(ARGS[5])
println(ARGS)
flush(stdout)



println("getting paths")
flush(stdout)
paths = valid_paths(dir, fname_pre)
res_paths = Tuple{String, String, Int64}[]
for i in 1:length(paths)
    path = paths[i]
    #TODO: config from input
    in_files = filter(x->startswith(x,fname_pre), readdir(path))
    for in_file in in_files
        push!(res_paths, (joinpath(path,in_file),joinpath(path,"../lDGA_julia/config.toml"),parse(Int,kIt)))
    end
end
println("constructed list of files: ")
map(p->println(p),res_paths)
flush(stdout)
res = pmap(x->E_pot_test(x...), res_paths)
println("------")
df = DataFrame(β = Float64[], U = Float64[], Nk = Int[], tc=String[],
               λsp = Float64[], λnew_sp = Float64[], λnew_ch = Float64[],
               EPot_direct_DMFT = Float64[], EPot_GS_DMFT = Float64[],
               EPot_chi_λ0 = Float64[], EPot_GS_λ0 = Float64[],
               EPot_chi_λsp = Float64[], EPot_GS_λsp = Float64[],
               EPot_chi_λspch = Float64[], EPot_GS_λspch = Float64[])
for (j,r) in enumerate(res)
    push!(df, r)
end
jldopen(out_fname,"a+") do f_out
    f_out["df"] = df
end
        # U = r[1]
        # beta = r[2]
        # Nk = r[3]
        # if !haskey(f_out, "$U@$beta")
        #     f_out["$U@$beta/beta"] = beta
        # f_out["$U@$beta/U"] = U
        # end
        # f_out["$U@$beta/$Nk/Nk"] = Nk
        # f_out["$U@$beta/$Nk/λsp"] = r[4]
        # f_out["$U@$beta/$Nk/λnew_sp"] = r[5]
        # f_out["$U@$beta/$Nk/λnew_ch"] = r[6]
        # f_out["$U@$beta/$Nk/E_pot_DMFT_chi"] = r[7]
        # f_out["$U@$beta/$Nk/E_pot_DMFT_sig"] = r[8]
        # f_out["$U@$beta/$Nk/E_pot_l0_chi"] = r[9]
        # f_out["$U@$beta/$Nk/E_pot_l0_sig"] = r[10]
        # f_out["$U@$beta/$Nk/E_pot_lsp_chi"] = r[11]
        # f_out["$U@$beta/$Nk/E_pot_lsp_sig"] = r[12]
        # f_out["$U@$beta/$Nk/E_pot_lspch_chi"] = r[13]
        # f_out["$U@$beta/$Nk/E_pot_lspch_sig"] = r[14]
    # end
