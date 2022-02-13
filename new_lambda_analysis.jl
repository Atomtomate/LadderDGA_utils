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
                        kG, mP::ModelParameters, sP::SimulationParameters)
    println("starting $(mP)")
    println(stderr,"starting $(mP)")
    flush(stdout)
    flush(stderr)
    # Prepare data
    res = Array{Float64, 2}(undef, length(λch_range), 6)
    ωindices = (sP.dbg_full_eom_omega) ? (1:size(nlQ_sp.χ,2)) : intersect(nlQ_sp.usable_ω, nlQ_ch.usable_ω)
    lur = length(ωindices)
    νmax = floor(Int,5*lur/16)
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
            E_kin, E_pot = calc_E(Σ_ladder.parent, kG, mP, sP)
            χupup_ω = subtract_tail(0.5 * kintegrate(kG,nlQ_ch_λ.χ .+ nlQ_sp_λ.χ,1)[1,ωindices], mP.Ekin_DMFT, iωn)
            χupdo_ω = 0.5 * kintegrate(kG,nlQ_ch_λ.χ .- nlQ_sp_λ.χ,1)[1,ωindices]
            lhs_c1 = real(sum(χupup_ω))/mP.β - mP.Ekin_DMFT*mP.β/12
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

function filter_usable_λsp_of_λch(λch_range, λsp_of_λch_data; max_λsp=Inf)
    #TODO: old version., why findmax??? 
    #tmp[isnan.(tmp)] .= 0.0
    #tmp[tmp .> max_λsp] .= 0.0
    #findmax(tmp)[2]:length(λch_range)
    tmp = deepcopy(λsp_of_λch_data)
    tmp = findall(x-> !isnan(x) && x < max_λsp, tmp)
end
