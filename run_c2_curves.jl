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
@everywhere Pkg.activate("/scratch/projects/hhp00048/codes/LadderDGA.jl")
@everywhere using LadderDGA
@everywhere using Dispersions
@everywhere using TOML
@everywhere using JLD2
@everywhere using Logging
@everywhere using FFTW
@everywhere using Setfield
@everywhere using NLsolve

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

        fine_grid, λnew_nlsolve = if hasfield(typeof(λ_int), :f_converged) && λ_int.f_converged
            println("U/β ($(mP.U)/$(mP.β)): nlsolve converged, found: $(λ_int.zero)")
            flush(stdout)
            range(λ_int.zero[2] - 0.1*λ_int.zero[2], λ_int.zero[2] + 0.1*λ_int.zero[2], 6), λ_int.zero
        else
            println("U/β ($(mP.U)/$(mP.β)): nlsolve not converged")
            flush(stdout)
            [], [NaN, NaN]
        end
        λch_range, spOfch = λsp_of_λch(in_f["nlQ_sp"], in_f["nlQ_ch"], kG, mP, sP, λsp_max=5.0, λch_max=500.0, n_λch=30-length(fine_grid), fine_grid=fine_grid)
        res = c2_along_λsp_of_λch(λch_range, spOfch, in_f["nlQ_sp"], in_f["nlQ_ch"], in_f["gLoc_fft"], λ₀, kG, mP, sP)


        nlQ_sp_λ = deepcopy(in_f["nlQ_sp"])
        nlQ_ch_λ = deepcopy(in_f["nlQ_ch"])
        ωindices = (1:size(nlQ_ch_λ.χ,2))
        νmax = minimum([sP.n_iν,floor(Int,3*length(ωindices)/8)])
        println("U/β ($(mP.U)/$(mP.β)): νmax = $νmax")
        # DMFT Energies
        Epot_DMFT_0 = mP.Epot_DMFT
        Σ_lDGA = calc_Σ(nlQ_sp_λ,nlQ_ch_λ,λ₀,in_f["gLoc_fft"],kG,mP, sP)
        Ekin_DMFT_1, Epot_DMFT_1 = calc_E(Σ_lDGA.parent[:,1:νmax], kG, mP, sP)
        Epot_DMFT_2 = real(0.5 * mP.U * sum(kintegrate(kG,nlQ_ch_λ.χ .- nlQ_sp_λ.χ,1)) / mP.β + mP.U * mP.n^2/4)
        println("U/β ($(mP.U)/$(mP.β)): E_pot DMFT: $(Epot_DMFT_1) / $(Epot_DMFT_2)")

        # Old λ Energies
        λsp_old = in_f["λsp_old"]
        nlQ_sp_λ.χ = LadderDGA.χ_λ(in_f["nlQ_sp"].χ, λsp_old)
        nlQ_sp_λ.λ = λsp_old
        Σ_lDGA_λsp = calc_Σ(nlQ_sp_λ, nlQ_ch_λ, λ₀, in_f["gLoc_fft"], kG, mP, sP)
        Ekin_lDGA_λsp_1, Epot_lDGA_λsp_1 = calc_E(Σ_lDGA_λsp.parent[:,1:νmax], kG, mP, sP)
        Epot_lDGA_λsp_2 = real(0.5 * mP.U * sum(kintegrate(kG,nlQ_ch_λ.χ .- nlQ_sp_λ.χ,1)) / mP.β + mP.U * mP.n^2/4)
        println("U/β ($(mP.U)/$(mP.β)): E_pot λsp: $(Epot_lDGA_λsp_1) / $(Epot_lDGA_λsp_2)")

        # New λ Energies
        λch_coarse = find_zero(res[:,2], res[:,5] .- res[:,6])
        LadderDGA.χ_λ!(nlQ_ch_λ.χ,in_f["nlQ_ch"].χ, λch_coarse)
        nlQ_ch_λ.λ = λch_coarse
        λsp_coarse = LadderDGA.λ_correction(:sp, in_f["imp_density"], in_f["nlQ_sp"], nlQ_ch_λ,
                                         in_f["gLoc_fft"], λ₀, kG, mP, sP)
        nls_res = LadderDGA.extended_λ(in_f["nlQ_sp"], in_f["nlQ_ch"],in_f["gLoc_fft"], λ₀, kG, mP, sP,
                    iterations=400, ftol=1e-9, x₀ = [λsp_coarse, λch_coarse])

        λsp_fine, λch_fine = hasfield(typeof(λ_int), :zero) ? nls_res.zero : nls_res
        λsp_new, λch_new = if sum(abs.([λsp_fine, λch_fine] .- [λsp_coarse, λch_coarse])) < 0.1
            λsp_fine, λch_fine
        else
            λsp_coarse, λch_coarse
        end
        LadderDGA.χ_λ!(nlQ_ch_λ.χ,in_f["nlQ_ch"].χ, λch_new)
        nlQ_ch_λ.λ = λch_new
        LadderDGA.χ_λ!(nlQ_sp_λ.χ,in_f["nlQ_sp"].χ, λsp_new)
        nlQ_sp_λ.λ = λsp_new
        Σ_lDGA_λspch = calc_Σ(nlQ_sp_λ, nlQ_ch_λ, λ₀, in_f["gLoc_fft"], kG, mP, sP)
        Ekin_lDGA_λspch_1, Epot_lDGA_λspch_1 = calc_E(Σ_lDGA_λspch.parent[:,1:νmax], kG, mP, sP)
        Epot_lDGA_λspch_2 = real(0.5 * mP.U * sum(kintegrate(kG,nlQ_ch_λ.χ .- nlQ_sp_λ.χ,1)) / mP.β + mP.U * mP.n^2/4)
        println("U/β ($(mP.U)/$(mP.β)): c2 done. found λsp, λch (coarse/fine) = $(λsp_coarse)/$(λsp_new), $(λch_coarse)/$(λch_new).\nEpot: $(Epot_lDGA_λspch_1) / $(Epot_lDGA_λspch_2)")
        flush(stdout)
        # E_pot
        return res, λsp_new, λch_new, λnew_nlsolve[1], λnew_nlsolve[2], Epot_DMFT_0, Epot_DMFT_1, Epot_DMFT_2, Epot_lDGA_λsp_1, Epot_lDGA_λsp_2, Epot_lDGA_λspch_1, Epot_lDGA_λspch_2, Ekin_DMFT_1, Ekin_lDGA_λsp_1, Ekin_lDGA_λspch_1
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
            f_out["$U@$beta/$Nk/Epot_DMFT_0"] = res[j][6]
            f_out["$U@$beta/$Nk/Epot_DMFT_1"] = res[j][7]
            f_out["$U@$beta/$Nk/Epot_DMFT_2"] = res[j][8]
            f_out["$U@$beta/$Nk/Epot_lDGA_λsp_1"] = res[j][9]
            f_out["$U@$beta/$Nk/Epot_lDGA_λsp_2"] = res[j][10]
            f_out["$U@$beta/$Nk/Epot_lDGA_λspch_1"] = res[j][11]
            f_out["$U@$beta/$Nk/Epot_lDGA_λspch_2"] = res[j][12]
            f_out["$U@$beta/$Nk/Ekin_DMFT_0"] = res[j][13]
            f_out["$U@$beta/$Nk/Ekin_DMFT_1"] = res[j][14]
            f_out["$U@$beta/$Nk/Ekin_DMFT_2"] = res[j][15]
        end
    end
end
