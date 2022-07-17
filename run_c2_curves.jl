"""
    Description:
Read through directory tree given by ARG1 and search for chi_lDGA output starting with ARG2. 
Calculate condition for lambda_ch as function of lambda_ch (root of this function is the lDGA 
value of lambda_ch). ARG4 is the number of k-points.
Results are written to directory ARG3.
use `collect_c2_curve_results.jl` to combine the results into a single file.

    Example:
julia run_c2_curves.jl /scratch/projects/hhp00048/lDGA/PD/data lDGA_ res 
"""

using Distributed
@everywhere using Pkg
@everywhere Pkg.activate("/scratch/projects/hhp00048/codes/LadderDGA.jl")
@everywhere using LadderDGA
@everywhere using Dispersions
@everywhere using TOML
@everywhere using JLD2, CodecZlib
@everywhere using Logging
@everywhere using FFTW
@everywhere using Setfield
@everywhere using NLsolve
using Base64

include("helpers/run_lDGA_dir.jl")
@everywhere include("new_lambda_analysis.jl")

@everywhere dir = $(ARGS[1])
@everywhere fname_pre = $(ARGS[2])
@everywhere out_dir = $(ARGS[3])
println(ARGS)
flush(stdout)
#TODO: ARGS[4] is a workaround, because some result files to not have the proper format yet

@everywhere function run_c2_curves(in::Tuple{String,String})
    in_file, out_file = in
    if isfile(out_file)
        println("file $out_file exists. Skipping computation.")
        return 1
    end
    out = IOBuffer()
    println(out,"starting with $in_file")
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
            println(out,"U/β ($(mP.U)/$(mP.β)): nlsolve converged, found: $(λ_int.zero)")
            range(λ_int.zero[2] - 0.1*λ_int.zero[2], λ_int.zero[2] + 0.1*λ_int.zero[2], 6), λ_int.zero
        else
            println(out,"U/β ($(mP.U)/$(mP.β)): nlsolve not converged")
            [], [NaN, NaN]
        end
        λch_range, spOfch = λsp_of_λch(in_f["nlQ_sp"], in_f["nlQ_ch"], kG, mP, sP, λsp_max=6.0, λch_max=mP.U*120.0, n_λch=30-length(fine_grid), fine_grid=fine_grid)
        println(out,"U/β ($(mP.U)/$(mP.β)): λch_range from = $(minimum(λch_range)) to $(maximum(λch_range)) with $(length(λch_range)) steps.\n
with λsp from = $(minimum(spOfch)) to $(maximum(spOfch)).")

        res = c2_along_λsp_of_λch(λch_range, spOfch, in_f["nlQ_sp"], in_f["nlQ_ch"], in_f["gLoc_fft"], λ₀, kG, mP, sP)


        nlQ_sp_λ = deepcopy(in_f["nlQ_sp"])
        nlQ_ch_λ = deepcopy(in_f["nlQ_ch"])
        ωindices = (1:size(nlQ_ch_λ.χ,2))
        νmax = minimum([sP.n_iν,floor(Int,3*length(ωindices)/8)])
        println(out,"U/β ($(mP.U)/$(mP.β)): νmax = $νmax")
        # DMFT Energies
        Epot_DMFT_0 = mP.Epot_DMFT
        Ekin_DMFT_0 = mP.Ekin_DMFT
        Σ_lDGA = calc_Σ(nlQ_sp_λ,nlQ_ch_λ,λ₀,in_f["gLoc_fft"],kG,mP, sP)
        Ekin_DMFT_1, Epot_DMFT_1 = calc_E(Σ_lDGA.parent[:,1:νmax], kG, mP, sP)
        Epot_DMFT_2 = real(0.5 * mP.U * sum(kintegrate(kG,nlQ_ch_λ.χ .- nlQ_sp_λ.χ,1)) / mP.β + mP.U * mP.n^2/4)
        println(out,"U/β ($(mP.U)/$(mP.β)): E_pot DMFT: $(Epot_DMFT_1) / $(Epot_DMFT_2)")

        # Old λ Energies
        λsp_old = in_f["λsp_old"]
        nlQ_sp_λ.χ = LadderDGA.χ_λ(in_f["nlQ_sp"].χ, λsp_old)
        nlQ_sp_λ.λ = λsp_old
        Σ_lDGA_λsp = calc_Σ(nlQ_sp_λ, nlQ_ch_λ, λ₀, in_f["gLoc_fft"], kG, mP, sP)
        Ekin_lDGA_λsp_1, Epot_lDGA_λsp_1 = calc_E(Σ_lDGA_λsp.parent[:,1:νmax], kG, mP, sP)
        Epot_lDGA_λsp_2 = real(0.5 * mP.U * sum(kintegrate(kG,nlQ_ch_λ.χ .- nlQ_sp_λ.χ,1)) / mP.β + mP.U * mP.n^2/4)
        println(out,"U/β ($(mP.U)/$(mP.β)): E_pot λsp: $(Epot_lDGA_λsp_1) / $(Epot_lDGA_λsp_2)")

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
        println(out,"U/β ($(mP.U)/$(mP.β)): c2 done. found λsp, λch (coarse/fine) = $(λsp_coarse)/$(λsp_new), $(λch_coarse)/$(λch_new).\nEpot: $(Epot_lDGA_λspch_1) / $(Epot_lDGA_λspch_2)")
        # E_pot
        println(out, "Writing to $out_file")
        log = String(take!(out))
        println(log)

        jldsave(out_file; res, λsp_new, λch_new, Epot_DMFT_0, Epot_DMFT_1, Epot_DMFT_2, Epot_lDGA_λsp_1, Epot_lDGA_λsp_2, Epot_lDGA_λspch_1, Epot_lDGA_λspch_2, Ekin_DMFT_0, Ekin_DMFT_1, Ekin_lDGA_λsp_1, Ekin_lDGA_λspch_1,log,in_file)
        return 0
    end
end

println("getting paths")
paths = valid_paths(dir, fname_pre)
res_paths = Tuple{String,String}[]

!isdir(out_dir) && mkdir(out_dir)

println("constructing list of files: ")
for i in 1:length(paths)
    path = paths[i]
    #TODO: config from input
    in_files = filter(x->startswith(x,fname_pre), readdir(path))
    for in_file in in_files
        inpath  = joinpath(path,in_file)
        outpath = joinpath(out_dir,base64encode(inpath)*".jld2")
        push!(res_paths, (inpath,outpath))
        println(inpath,"\n",outpath,"\n============")
    end
end
res = nprocs() > 1 ? pmap(run_c2_curves, res_paths; on_error=ex->nothing) : map(run_c2_curves, res_paths)
println("====================\nDONE\n====================\nSTATUS:")
println(res)
