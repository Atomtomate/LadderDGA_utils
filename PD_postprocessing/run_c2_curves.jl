"""
    Description:
Read through directory tree given by ARG1 and search for chi_lDGA output starting with ARG2. 
Calculate condition for lambda_ch as function of lambda_ch (root of this function is the lDGA 
value of lambda_ch).
Results are written to ARG3.

    Example:
julia run_c2_curves.jl /scratch/projects/hhp00048/lDGA/PD/data lDGA_ res
"""

using Distributed
@everywhere using Pkg
@everywhere Pkg.activate("/scratch/projects/hhp00048/codes/LadderDGA.jl",io=devnull)
@everywhere using LadderDGA
@everywhere using TOML
@everywhere using JLD2
@everywhere using Logging

include("helpers/run_lDGA_dir.jl")
@everywhere include("/scratch/projects/hhp00048/codes/LadderDGA_utils/new_lambda_analysis.jl")

@everywhere dir = $(ARGS[1])
@everywhere fname_pre = $(ARGS[2])
@everywhere out_fname = $(ARGS[3])
@everywhere Nk = $(ARGS[4])
println(ARGS)
flush(stdout)

#if !isdir(out_fname*"_logs")
#    mkdir(out_fname*"_logs")
#end
#out_fname*"logs/w$(myid())_.log","w"
#for p in workers()
#    @spawnat p LadderDGA.redirect_logs(out_fname*"logs/w$(myid())_.log")
#end
#TODO: ARGS[4] is a workaround, because some result files to not have the proper format yet

@everywhere function run_c2_curves(in_file)
    Nk = 10
    println("starting with $in_file")
    flush(stdout)
    jldopen(in_file, "r") do in_f
        cfg = in_f["config"]
        mP, sP, env, kGridsStr = readConfig(cfg);
        println("config file read of: $in_file")
        flush(stdout)
        #TODO: cfg from in_file, kIteration also
        impQ_sp, impQ_ch, _, _, kG, _, gLoc_fft, Σ_loc, FUpDo = setup_LDGA(("3Dsc-0.2041241452319315", Nk), mP, sP, env);
        λch_range, spOfch = λsp_of_λch(in_f["nlQ_sp"], in_f["nlQ_ch"], kG, mP, sP, max_λsp=30.0, λch_max=25.0, n_λch=500)
        res = c2_along_λsp_of_λch(λch_range, spOfch, in_f["nlQ_sp"], in_f["nlQ_ch"],
                            in_f["bubble"], in_f["Sigma_loc"], Σ_loc,
                            gLoc_fft, FUpDo, kG, mP, sP)
        λch = find_zero(res[:,2], res[:,5] .- res[:,6])
        nlQ_ch_λ = deepcopy(in_f["nlQ_ch"])
        nlQ_ch_λ.χ = LadderDGA.χ_λ(in_f["nlQ_ch"].χ, λch)
        nlQ_ch_λ.λ = λch
        λsp_new = LadderDGA.λ_correction(:sp,impQ_sp, impQ_ch, FUpDo, Σ_loc, in_f["Sigma_loc"], in_f["nlQ_sp"], nlQ_ch_λ,
                               in_f["bubble"], gLoc_fft, kG, mP, sP)
        λ_new = LadderDGA.λ_correction(:sp_ch,impQ_sp, impQ_ch, FUpDo, Σ_loc, in_f["Sigma_loc"], in_f["nlQ_sp"], in_f["nlQ_ch"],
                               in_f["bubble"], gLoc_fft, kG, mP, sP)
        println("done found lsp,lch = $(λsp_new), $(λch)")
        flush(stdout)
        return mP.U,mP.β,Nk,res, λsp_new, λch, λ_new
    end
end

println("getting paths")
flush(stdout)
paths = valid_paths(dir, fname_pre)
res_paths = String[]
for i in 1:length(paths)
    path = paths[i]
    #TODO: config from input
    in_files = filter(x->startswith(x,fname_pre), readdir(path))
    for in_file in in_files
        println("combining $(joinpath(path,in_file))")
        flush(stdout)
        push!(res_paths, joinpath(path,in_file))
    end
end
println("constructed list of files: ")
println(res_paths)
flush(stdout)
#println(res_f)
#res = fetch.(res_f)
#push!(res_f, @spawnat :any run_c2_curves(joinpath(path,in_file)))
res = pmap(run_c2_curves, res_paths)
println("------")
#println(res)


jldopen(out_fname,"a+") do f_out
    for (j,p) in enumerate(res_paths)
        jldopen(p, "r") do f_in
            #println(keys(f_in))
            U = res[j][1]
            beta = res[j][2]
            Nk = res[j][3]
            if !haskey(f_out, "$U@$beta")
                f_out["$U@$beta/beta"] = beta
            f_out["$U@$beta/U"] = U
            end
            f_out["$U@$beta/$Nk/Nk"] = Nk
            f_out["$U@$beta/$Nk/res"] = res[j][4]
            f_out["$U@$beta/$Nk/λsp"] = res[j][5]
            f_out["$U@$beta/$Nk/λch"] = res[j][6]
            f_out["$U@$beta/$Nk/λnew"] = res[j][7]
            f_out["$U@$beta/$Nk/λsp_old"] = f_in["λsp_old"]
            f_out["$U@$beta/$Nk/χsp"] = f_in["nlQ_sp"].χ
            f_out["$U@$beta/$Nk/χch"] = f_in["nlQ_ch"].χ
        end
    end
end
