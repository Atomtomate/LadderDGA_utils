using Dispersions

function filter_usable_λsp_of_λch(λch_range, inp; max_λsp=Inf)
    tmp = deepcopy(inp)
    tmp[isnan.(tmp)] .= 0.0
    tmp[tmp .> max_λsp] .= 0.0
    findmax(tmp)[2]:length(λch_range)
end



function E_pot_test(bubble, impQ_sp, impQ_ch, nlQ_sp, nlQ_ch, Σ_ladder_loc, kn, gLoc, gLoc_fft, Σ_loc, FUpDo, mP, sP)
    res = try
    #TODO: DO NOT HARDCODE LATTICE TYPE!!!
    kG = gen_kGrid("3Dsc-0.2041241452319315", kn)
    mP.Epot_DMFT
    gL = transpose(LadderDGA.flatten_2D(gLoc[sP.fft_offset:(sP.fft_offset+299)]));
    sL = Array{ComplexF64,2}(undef, length(kG.kMult), 300)
    for i in 1:300
        sL[:,i] .= Σ_loc[i]
    end

    Ek1, Ep1 = calc_E(sL, kG, mP, sP)

    # local
    #bubbleLoc = calc_bubble(gImp, kGridLoc, mP, sP);
    #locQ_sp = calc_χ_trilex(impQ_sp.Γ, bubbleLoc, kGridLoc, mP.U, mP, sP);
    #locQ_ch = calc_χ_trilex(impQ_ch.Γ, bubbleLoc, kGridLoc, -mP.U, mP, sP);
    #Σ_ladder_loc = calc_Σ(locQ_sp, locQ_ch, bubbleLoc, gImp, FUpDo, kGridLoc, mP, sP)

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
    ωindices = intersect(nlQ_sp.usable_ω, nlQ_ch.usable_ω)
    νmax = floor(Int,length(ωindices)/3)

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
    χupdo_ω = real.(kintegrate(kG, nlQ_ch.χ .- nlQ_sp.χ, 1)[1,ωindices])
    E_pot3 = real(sum(χupdo_ω)/(2*mP.β)) + (mP.n/2)^2

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
