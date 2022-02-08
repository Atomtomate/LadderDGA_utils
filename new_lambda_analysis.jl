using Dispersions
using LadderDGA

function find_zero(λch_vals::AbstractVector, c2_curve::AbstractVector)
    yvals = c2_curve
    xvals = λch_vals
    sc_ind = findfirst(x->x<0,sign.(yvals))
    (sc_ind == nothing) && return Inf
    (sc_ind == 1) && return -Inf
    y1 = yvals[sc_ind-1]
    y2 = yvals[sc_ind]
    x1 = xvals[sc_ind-1]
    x2 = xvals[sc_ind]
    m = (y2-y1)/(x2-x1)
    x0 = x1 - y1/m
end

function find_epot(λch_vals, c2_curve, res)
    rrr = []
    yvals = c2_curve
    xvals = λch_vals
    sc_ind = findfirst(x->x<0,sign.(yvals))
    (sc_ind == nothing) && return [NaN, NaN, NaN, NaN]
    (sc_ind == 1) && return [NaN, NaN, NaN, NaN]
    return res[sc_ind-1,5], res[sc_ind-1,6], res[sc_ind,5], res[sc_ind,6]
end

function λsp_of_λch(nlQ_sp::NonLocalQuantities, nlQ_ch::NonLocalQuantities, kG, mP, sP; max_λsp=10.0, λch_max=10.0, n_λch=50, fine_grid=[])
    ωindices = (sP.dbg_full_eom_omega) ? (1:size(nlQ_sp.χ,2)) : intersect(nlQ_sp.usable_ω, nlQ_ch.usable_ω)
    iωn = 1im .* 2 .* (-sP.n_iω:sP.n_iω)[ωindices] .* π ./ mP.β
    nh  = ceil(Int64, size(nlQ_sp.χ,2)/2)
    χsp_min    = -minimum(1 ./ real.(nlQ_sp.χ[:,nh]))
    χch_min    = -minimum(1 ./ real.(nlQ_ch.χ[:,nh]))

    λch_range = Float64.(sort(union(range(χch_min,λch_max,n_λch), [0], fine_grid)))
    spOfch_max_nl = zeros(length(λch_range))

    χsp_nλ_r = real.(deepcopy(nlQ_sp.χ[:,ωindices]))
    χch_nλ_r = real.(deepcopy(nlQ_ch.χ[:,ωindices]))
    χch_λ = similar(χch_nλ_r)

    for (λi,λchi) in enumerate(λch_range)
        χ_λ!(χch_λ, χch_nλ_r, λchi)
        χch_ω = kintegrate(kG, χch_λ, 1)[1,:]
        χch_ω_sub = subtract_tail(χch_ω, mP.Ekin_DMFT, iωn)
        χch_sum = real(sum(χch_ω_sub))/mP.β - mP.Ekin_DMFT*mP.β/12
        rhs_val = (mP.n) * (1 - mP.n/2) - χch_sum
        λ = λsp(χsp_nλ_r, iωn, mP.Ekin_DMFT, rhs_val, kG, mP)
        spOfch_max_nl[λi] = λ
    end;

    λch_range_filtered = filter_usable_λsp_of_λch(λch_range, spOfch_max_nl; max_λsp=max_λsp)
    λch_range_f = λch_range[λch_range_filtered]
    spOfch_f = spOfch_max_nl[λch_range_filtered]
    return λch_range_f, spOfch_f
end

function c2_along_λsp_of_λch(λch_range::AbstractArray{Float64,1}, spOfch::AbstractArray{Float64,1},
                        nlQ_sp::NonLocalQuantities, nlQ_ch::NonLocalQuantities,
                        Gνω::GνqT, λ₀::AbstractArray{ComplexF64,3}, 
                        kG::ReducedKGrid, mP::ModelParameters, sP::SimulationParameters)
    println("starting $(mP)")
    println(stderr,"starting $(mP)")
    flush(stdout)
    flush(stderr)
    # Prepare data
    res = Array{Float64, 2}(undef, length(λch_range), 6)
    ωindices = (sP.dbg_full_eom_omega) ? (1:size(nlQ_sp.χ,2)) : intersect(nlQ_sp.usable_ω, nlQ_ch.usable_ω)
    lur = length(ωindices)
    νmax = floor(Int,3*lur/8)
    if νmax > 5
        iωn = 1im .* 2 .* (-sP.n_iω:sP.n_iω)[ωindices] .* π ./ mP.β
        nlQ_sp_λ = deepcopy(nlQ_sp)
        nlQ_ch_λ = deepcopy(nlQ_ch)
        for i in 1:length(λch_range)
            λ = [spOfch[i], λch_range[i]]
            χ_λ!(nlQ_sp_λ.χ, nlQ_sp.χ, λ[1])
            χ_λ!(nlQ_ch_λ.χ, nlQ_ch.χ, λ[2])
            Σ_ladder = calc_Σ(nlQ_sp_λ, nlQ_ch_λ, λ₀, Gνω, kG, mP, sP)[:,0:νmax]
            #Σ_ladder = Σ_loc_correction(Σ_ladder, Σ_ladderLoc, Σ_loc);
            
            χupup_ω = subtract_tail(0.5 * kintegrate(kG,nlQ_ch_λ.χ .+ nlQ_sp_λ.χ,1)[1,ωindices], mP.Ekin_DMFT, iωn)
            χupdo_ω = 0.5 * kintegrate(kG,nlQ_ch_λ.χ .- nlQ_sp_λ.χ,1)[1,ωindices]
            E_kin, E_pot = calc_E(Σ_ladder.parent, kG, mP, sP)
            lhs_c1 = real(sum(χupup_ω))/mP.β - E_kin*mP.β/12
            lhs_c2 = real(sum(χupdo_ω))/mP.β
            rhs_c1 = mP.n/2 * (1 - mP.n/2)
            rhs_c2 = E_pot/mP.U - (mP.n/2) * (mP.n/2)
            res[i,:] = [λ[1] λ[2] lhs_c1 rhs_c1 lhs_c2 rhs_c2]
        end
    end
    return res
end

function new_λ_from_c2(c2_res, imp_dens, nlQ_sp, nlQ_ch, locQ_sp, gLoc_fft, λ₀, kG, mP, sP)
    λsp, λch = if size(c2_res,1) >= 1
        λch = find_zero(c2_res[:,2], c2_res[:,5] .- c2_res[:,6])
        nlQ_ch_λ = deepcopy(nlQ_ch)
        nlQ_ch_λ.χ = LadderDGA.χ_λ(nlQ_ch.χ, λch)
        nlQ_ch_λ.λ = λch
        λsp = LadderDGA.λ_correction(:sp, imp_dens, nlQ_sp, nlQ_ch_λ, gLoc_fft, λ₀, kG, mP, sP)
        λsp, λch
    else
        NaN, NaN
    end
    λsp, λch
end

function filter_usable_λsp_of_λch(λch_range, inp; max_λsp=Inf)
    tmp = deepcopy(inp)
    tmp[isnan.(tmp)] .= 0.0
    tmp[tmp .> max_λsp] .= 0.0
    findmax(tmp)[2]:length(λch_range)
end



function E_pot_test(λsp, λnew, bubble, nlQ_sp, nlQ_ch, Σ_ladder_loc, kn, gLoc, gLoc_fft, Σ_loc, Fsp, kG, mP, sP)
    res = try
    #TODO: DO NOT HARDCODE LATTICE TYPE!!!
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
    #Σ_ladder_loc = calc_Σ(locQ_sp, locQ_ch, bubbleLoc, gImp, Fsp, kGridLoc, mP, sP)
    # lambda
    nlQ_sp_λsp = χ_λ(nlQ_sp, λsp)
    
    nlQ_sp_λnew = χ_λ(nlQ_sp, λnew[1])
    nlQ_ch_λnew = χ_λ(nlQ_ch, λnew[2])
    
    # Self energies
    ωindices = intersect(nlQ_sp.usable_ω, nlQ_ch.usable_ω)
    νmax = floor(Int,length(ωindices)/3)

    Σ_ladder = calc_Σ(nlQ_sp, nlQ_ch, bubble, gLoc_fft, Fsp, kG, mP, sP)[:,1:νmax]
    Σ_ladder = Σ_loc_correction(Σ_ladder, Σ_ladder_loc, Σ_loc);
    E_kin_lDGA,E_pot_lDGA = LadderDGA.calc_E(Σ_ladder, kG, mP, sP)

    Σ_ladder_λsp = calc_Σ(nlQ_sp_λsp, nlQ_ch, bubble, gLoc_fft, Fsp, kG, mP, sP)[:,1:νmax]
    Σ_ladder_λsp = Σ_loc_correction(Σ_ladder_λsp, Σ_ladder_loc, Σ_loc);
    E_kin_lDGA_λsp, E_pot_lDGA_λsp = LadderDGA.calc_E(Σ_ladder_λsp, kG, mP, sP)

    Σ_ladder_λnew = calc_Σ(nlQ_sp_λnew, nlQ_ch_λnew, bubble, gLoc_fft, Fsp, kG, mP, sP)[:,1:νmax]
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
         mP.Epot_DMFT, Ep1, 
         E_pot_lDGA/mP.U, E_pot3, E_pot_lDGA_λsp/mP.U, E_pot3_λsp, E_pot_lDGA_λnew/mP.U, E_pot3_λnew,
         nlQ_sp.χ, nlQ_ch.χ]
    catch e
        println(stderr, "caught error: ", e)
        [mP.U, mP.β, kG.Ns, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN,
         Matrix{ComplexF64}[], Matrix{ComplexF64}[]]
    end
    return res
end
