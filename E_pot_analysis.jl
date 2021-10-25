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
@everywhere include("new_lambda_analysis.jl")

@everywhere dir = $(ARGS[1])
@everywhere fname_pre = $(ARGS[2])
@everywhere out_fname = $(ARGS[3])
@everywhere kIt = $(ARGS[4])
@everywhere tcType = $(ARGS[5])
println(ARGS)
flush(stdout)

@everywhere function E_pot_test(file, config, kIt)
    mP, sP, env, kGridsStr = readConfig(config)
    impQ_sp, impQ_ch, gImp, kGridLoc, kG, gLoc, gLoc_fft, Σ_loc, FUpDo = setup_LDGA(kGridsStr[kIt], mP, sP, env);
    res = try
    bubble, nlQ_sp, nlQ_ch = jldopen(file,"r") do f
        f["bubble"], f["nlQ_sp"], f["nlQ_ch"]
    end

    mP.Epot_DMFT
    gL = transpose(LadderDGA.flatten_2D(gLoc[sP.fft_offset:(sP.fft_offset+299)]));
    sL = Array{ComplexF64,2}(undef, length(kG.kMult), 300)
    for i in 1:300
        sL[:,i] .= Σ_loc[i]
    end

    Ek1, Ep1 = calc_E(sL, kG, mP, sP)
    ωindices = intersect(impQ_ch.usable_ω, impQ_ch.usable_ω)
    νmax = floor(Int,length(ωindices)/3)
    E_pot2 = real(sum((impQ_ch.χ_ω .- impQ_sp.χ_ω)[ωindices])/(2*mP.β))

    # local
    bubbleLoc = calc_bubble(gImp, kGridLoc, mP, sP);
    locQ_sp = calc_χ_trilex(impQ_sp.Γ, bubbleLoc, kGridLoc, mP.U, mP, sP);
    locQ_ch = calc_χ_trilex(impQ_ch.Γ, bubbleLoc, kGridLoc, -mP.U, mP, sP);
    Σ_ladder_loc = calc_Σ(locQ_sp, locQ_ch, bubbleLoc, gImp, FUpDo, kGridLoc, mP, sP)

    # lambda
    λsp = λ_correction(:sp, impQ_sp, impQ_ch, FUpDo, Σ_loc, Σ_ladder_loc, nlQ_sp,nlQ_ch,
                    bubble, gLoc_fft, kG, mP, sP)
    nlQ_sp_λsp = χ_λ(nlQ_sp, λsp)

    λnew = λ_correction(:sp_ch, impQ_sp, impQ_ch, FUpDo, Σ_loc, Σ_ladder_loc, nlQ_sp, nlQ_ch,
                    bubble, gLoc_fft, kG, mP, sP)
    λnew = 
        if λnew.f_converged 
            @info λnew;
            λnew.zero
        else 
            [NaN, NaN] 
        end
        
    #λnew = [NaN, NaN]
    
    nlQ_sp_λnew = χ_λ(nlQ_sp, λnew[1])
    nlQ_ch_λnew = χ_λ(nlQ_ch, λnew[2])
    
    # Self energies
    Σ_ladder = calc_Σ(nlQ_sp, nlQ_ch, bubble, gLoc_fft, FUpDo, kG, mP, sP)[:,1:νmax]
    Σ_ladder = Σ_loc_correction(Σ_ladder, Σ_ladder_loc, Σ_loc);
    E_kin_lDGA,E_pot_lDGA = LadderDGA.calc_E(Σ_ladder, kG, mP, sP)

    Σ_ladder_λsp = calc_Σ(nlQ_sp_λsp, nlQ_ch, bubble, gLoc_fft, FUpDo, kG, mP, sP)[:,1:νmax]
    Σ_ladder_λsp = Σ_loc_correction(Σ_ladder_λsp, Σ_ladder_loc, Σ_loc);
    E_kin_lDGA_λsp, E_pot_lDGA_λsp = LadderDGA.calc_E(Σ_ladder_λsp, kG, mP, sP)

    Σ_ladder_λnew = calc_Σ(nlQ_sp_λnew, nlQ_ch_λnew, bubble, gLoc_fft, FUpDo, kG, mP, sP)[:,1:νmax]
    Σ_ladder_λnew = Σ_loc_correction(Σ_ladder_λnew, Σ_ladder_loc, Σ_loc);
    E_kin_lDGA_λnew, E_pot_lDGA_λnew = LadderDGA.calc_E(Σ_ladder_λnew, kG, mP, sP)

    
    # chi
    ωindices = intersect(nlQ_sp.usable_ω, nlQ_ch.usable_ω)
    χupdo_ω = real.(kintegrate(kG, nlQ_ch.χ .- nlQ_sp.χ, 1)[1,ωindices])
    E_pot3 = real(sum(χupdo_ω)/(2*mP.β)) + (mP.n/2)^2

    ωindices = intersect(nlQ_sp.usable_ω, nlQ_ch.usable_ω)
    χupdo_λsp_ω = real.(kintegrate(kG, nlQ_ch.χ .- nlQ_sp_λsp.χ, 1)[1,ωindices])
    E_pot3_λsp = real(sum(χupdo_λsp_ω)/(2*mP.β)) + (mP.n/2)^2
    
    χupdo_λnew_ω = real.(kintegrate(kG, nlQ_ch_λnew.χ .- nlQ_sp_λnew.χ, 1)[1,ωindices])
    E_pot3_λnew = real(sum(χupdo_λnew_ω)/(2*mP.β)) + (mP.n/2)^2

    println("\n ============================================= \n")
    println(mP)
    println("U<n_up n_do >                 = $(mP.Epot_DMFT)")
    println("∑_ν G_loc Σ_loc               = $(Ep1)")
    println("1/(β U) ∑_ν G_{λ=0} Σ_{λ=0}   = $(E_pot_lDGA/mP.U)")
    println("∑_{ωq} χ^{λ=0}_{updo} + n²/4  = $(E_pot3)")
    #println("∑_{ωq} χ^imp_{updo}   = $(E_pot2)")
    println("1/(β U) ∑_ν G_λsp Σ_λsp       = $(E_pot_lDGA_λsp/mP.U)")
    println("∑_{ωq} χ^{λsp}_{updo} + n²/4  = $(E_pot3_λsp)")
    println("1/(β U) ∑_ν G_λnew Σ_λnew     = $(E_pot_lDGA_λnew/mP.U)")
    println("∑_{ωq} χ^{λnew}_{updo} + n²/4 = $(E_pot3_λnew)")
    println(" ============================================= \n")
        [mP.U, mP.β, kG.Ns, 
         λsp, λnew[1], λnew[2], mP.Epot_DMFT, Ep1, 
         E_pot_lDGA/mP.U, E_pot3, E_pot_lDGA_λsp/mP.U, E_pot3_λsp, E_pot_lDGA_λnew/mP.U, E_pot3_λnew,
         nlQ_sp.χ, nlQ_ch.χ]
    catch e
        println(stderr, "caught error: ", e)
        [mP.U, mP.β, kG.Ns, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN,
         Matrix{ComplexF64}[], Matrix{ComplexF64}[]]
    end
    return res
end


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
