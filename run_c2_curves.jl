"""
    Description:
Read through directory tree given by ARG1 and search for chi_lDGA output starting with ARG2. 
Calculate condition for lambda_ch as function of lambda_ch (root of this function is the lDGA 
value of lambda_ch). ARG4 is the number of k-points.
Results are written to ARG3.

    Example:
julia run_c2_curves.jl /scratch/projects/hhp00048/lDGA/PD/data lDGA_ res 
"""

using Distributed
@everywhere using Pkg
@everywhere Pkg.activate("/scratch/projects/hhp00048/codes/LadderDGA.jl",io=devnull)
@everywhere using LadderDGA
@everywhere using Dispersions
@everywhere using TOML
@everywhere using JLD2
@everywhere using Logging
@everywhere using FFTW
@everywhere using Setfield

include("helpers/run_lDGA_dir.jl")
@everywhere include("new_lambda_analysis.jl")

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
    println("starting with $in_file")
    flush(stdout)
    jldopen(in_file, "r") do in_f
        mP = in_f["mP"]
        sP = in_f["sP"]
        kG = in_f["kG"]
        p = plan_fft!(Array{ComplexF64}(undef, gridshape(kG)...), flags=FFTW.ESTIMATE, timelimit=Inf)
        kG = @set kG.fftw_plan = p
        λ₀ = in_f["λ₀_sp"]

        #TODO: or read this...
        λ_int = LadderDGA.extended_λ(in_f["nlQ_sp"], in_f["nlQ_ch"], in_f["gLoc_fft"], λ₀, kG, mP, sP, iterations=20, ftol=1e-5)

        fine_grid, λnew_nlsolve = if λ_int.f_converged
            println("nlsolve converged, found: $(λ_int.zero)")
            flush(stdout)
            range(λ_int.zero[2] - 0.1*λ_int.zero[2], λ_int.zero[2] + 0.1*λ_int.zero[2], 20), λ_int.zero
        else
            println("nlsolve not converged")
            flush(stdout)
            [], [NaN, NaN]
        end
        λch_range, spOfch = λsp_of_λch(in_f["nlQ_sp"], in_f["nlQ_ch"], kG, mP, sP, max_λsp=5.0, λch_max=30.0, n_λch=70-length(fine_grid), fine_grid=fine_grid)
        res = c2_along_λsp_of_λch(λch_range, spOfch, in_f["nlQ_sp"], in_f["nlQ_ch"], in_f["gLoc_fft"], λ₀, kG, mP, sP)
        λch = find_zero(res[:,2], res[:,5] .- res[:,6])
        nlQ_ch_λ = deepcopy(in_f["nlQ_ch"])
        nlQ_ch_λ.χ = LadderDGA.χ_λ(in_f["nlQ_ch"].χ, λch)
        nlQ_ch_λ.λ = λch
        λsp_new = LadderDGA.λ_correction(:sp, in_f["imp_density"], in_f["nlQ_sp"], nlQ_ch_λ,
                                         in_f["gLoc_fft"], λ₀, kG, mP, sP)
        println("c2 for β=$(mP.β), U=$(mP.U) done. found λsp, λch = $(λsp_new), $(λch)")
        flush(stdout)
        return res, λsp_new, λch, λnew_nlsolve[1], λnew_nlsolve[2]
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
res = nprocs() > 1 ? pmap(run_c2_curves, res_paths) : map(run_c2_curves, res_paths)
println("------")
#println(res)


jldopen(out_fname,"a+") do f_out
    for (j,p) in enumerate(res_paths)
        jldopen(p, "r") do f_in
            #println(keys(f_in))
            U = f_in["mP"].U
            beta = f_in["mP"].β
            Nk = if haskey(f_in, "Nk")
                f_in["Nk"]
            else
                parse(Int,ARGS[4])^3
            end
            if !haskey(f_out, "$U@$beta")
                f_out["$U@$beta/beta"] = beta
            f_out["$U@$beta/U"] = U
            end
            #res_j = run_c2_curves(p)
            #println(res)
            f_out["$U@$beta/$Nk/Nk"] = Nk
            f_out["$U@$beta/$Nk/res"] = res[j][1]
            f_out["$U@$beta/$Nk/λsp"] = res[j][2]
            f_out["$U@$beta/$Nk/λch"] = res[j][3]
            f_out["$U@$beta/$Nk/λsp_nlsolve"] = res[j][4]
            f_out["$U@$beta/$Nk/λch_nlsolve"] = res[j][5]
            f_out["$U@$beta/$Nk/χ_sp"] = f_in["nlQ_sp"].χ
            f_out["$U@$beta/$Nk/χ_ch"] = f_in["nlQ_ch"].χ
            f_out["$U@$beta/$Nk/lDGAlog"] = f_in["log"]
            f_out["$U@$beta/$Nk/χ_sp_usable"] = f_in["nlQ_sp"].usable_ω
            f_out["$U@$beta/$Nk/χ_ch_usable"] = f_in["nlQ_ch"].usable_ω
            f_out["$U@$beta/$Nk/λsp_old"] = f_in["λsp_old"]
        end
    end
end
