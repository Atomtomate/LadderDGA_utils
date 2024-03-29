function find_zero(λch_vals::AbstractVector, c2_curve::AbstractVector)
    yvals = c2_curve
    xvals = λch_vals
    sc_ind = findlast(x -> x == -2, diff(sign.(c2_curve)))
    (sc_ind == nothing) && return Inf
    (sc_ind == 1) && return -Inf
    y1 = yvals[sc_ind]
    y2 = yvals[sc_ind+1]
    x1 = xvals[sc_ind]
    x2 = xvals[sc_ind+1]
    m = (y2-y1)/(x2-x1)
    x0 = x1 - y1/m
end

function find_epot(λch_vals, c2_curve, res)
    rrr = []
    yvals = c2_curve
    xvals = λch_vals
    sc_ind = findlast(x -> x == -2, diff(sign.(c2_curve)))
    (sc_ind == nothing) && return [NaN, NaN, NaN, NaN]
    (sc_ind == 1) && return [NaN, NaN, NaN, NaN]
    return res[sc_ind,5], res[sc_ind,6], res[sc_ind+1,5], res[sc_ind+1,6]
end

function λsp_of_λch(χsp, χch, kG, mP, sP; λsp_max=10.0, λch_max=1000.0, n_λch=50, fine_grid=[],
                    range_spacing_exp=4.0)
    ωindices = (sP.dbg_full_eom_omega) ? (1:size(χsp.data,2)) : intersect(χsp.usable_ω, χch.usable_ω)
    iωn = 1im .* 2 .* (-sP.n_iω:sP.n_iω)[ωindices] .* π ./ mP.β
    nh  = ceil(Int64, size(χsp.data,2)/2)
    #TODO: find reason for extremely large χch_min at U>3
    χsp_min    = -minimum(1 ./ real.(χsp.data[:,nh]))
    χch_min    = -minimum(1 ./ real.(χch.data[:,nh]))
    χch_min = if χch_min > 500
        println("WARNING: found χ inv min = $χch_min. Resetting to -1!")
        -5.0
    else
        χch_min
    end
    @warn "found χsp_min: " χsp_min ", χch_min" χch_min

    # Construct log spacing (function is smooth for large λch
    range_shift = -χch_min
    r = range(χch_min+range_shift,(λch_max+range_shift)^(1/range_spacing_exp),n_λch).^range_spacing_exp .- range_shift
    λch_range = Float64.(sort(union(r, [0], fine_grid)))
    @warn "constructed range: " λch_range
    spOfch_max_nl = zeros(length(λch_range))

    χsp_nλ_r = real.(deepcopy(χsp.data[:,ωindices]))
    χch_nλ_r = real.(deepcopy(χch.data[:,ωindices]))
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
    @warn "lsp_of_ch" spOfch_max_nl
    #DBG: deactivated filter by setting min to -Inf
    λch_range_filtered = filter_usable_λsp_of_λch(λch_range, spOfch_max_nl, -Inf, -Inf; λsp_max=λsp_max, λch_max=λch_max)
    #λch_range_filtered = filter_usable_λsp_of_λch(λch_range, spOfch_max_nl, χsp_min, χch_min; λsp_max=λsp_max, λch_max=λch_max)
    λch_range_f = λch_range[λch_range_filtered]

    @warn "filtered range: " λch_range_f
    spOfch_f = spOfch_max_nl[λch_range_filtered]
    return λch_range_f, spOfch_f
end

function c2_along_λsp_of_λch(λch_range::AbstractArray{Float64,1}, spOfch::AbstractArray{Float64,1},
                        χsp::χT, χch::χT, γsp::γT, γch::γT,
                        Gνω::GνqT, λ₀::AbstractArray{ComplexF64,3}, 
                        kG, mP::ModelParameters, sP::SimulationParameters)
    println("starting $(mP)")
    println(stderr,"starting $(mP)")
    flush(stdout)
    flush(stderr)
    # Prepare data
    res = Array{Float64, 2}(undef, length(λch_range), 6)
    ωindices = (sP.dbg_full_eom_omega) ? (1:size(χsp.data,2)) : intersect(χsp.usable_ω, χch.usable_ω)
    lur = length(ωindices)
    νmax = floor(Int,5*lur/16)
    if νmax > 5
        iωn = 1im .* 2 .* (-sP.n_iω:sP.n_iω)[ωindices] .* π ./ mP.β
        χsp_λ = deepcopy(χsp)
        χch_λ = deepcopy(χch)
        for i in 1:length(λch_range)
            λ = [spOfch[i], λch_range[i]]
            χ_λ!(χsp_λ.χ, χsp.χ, λ[1])
            χ_λ!(χch_λ.χ, χch.χ, λ[2])
            Σ_ladder = calc_Σ(χsp, χch, γsp, γch, λ₀, Gνω, kG, mP, sP)[:,0:νmax]
            #Σ_ladder = Σ_loc_correction(Σ_ladder, Σ_ladderLoc, Σ_loc);
            E_kin, E_pot = calc_E(Σ_ladder.parent, kG, mP, sP)
            χupup_ω = subtract_tail(0.5 * kintegrate(kG,χch_λ.data .+ χsp_λ.data,1)[1,ωindices], mP.Ekin_DMFT, iωn)
            χupdo_ω = 0.5 * kintegrate(kG,χch_λ.data .- χsp_λ.data,1)[1,ωindices]
            lhs_c1 = real(sum(χupup_ω))/mP.β - mP.Ekin_DMFT*mP.β/12
            lhs_c2 = real(sum(χupdo_ω))/mP.β
            rhs_c1 = mP.n/2 * (1 - mP.n/2)
            rhs_c2 = E_pot/mP.U - (mP.n/2) * (mP.n/2)
            res[i,:] = [λ[1] λ[2] lhs_c1 rhs_c1 lhs_c2 rhs_c2]
        end
    end
    return res
end

function filter_usable_λsp_of_λch(λch_range, λsp_of_λch_data, χsp_min, χch_min; λsp_max=Inf, λch_max=Inf)
    #TODO: old version., why findmax??? 
    #tmp[isnan.(tmp)] .= 0.0
    #tmp[tmp .> max_λsp] .= 0.0
    #findmax(tmp)[2]:length(λch_range)
    χch_filter_indices = findall(x -> !isnan(x) && x < λch_max && x > χch_min, λch_range)
    χsp_filter_indices = findall(x -> !isnan(x) && x < λsp_max && x > χsp_min, λsp_of_λch_data)
    return sort(intersect(χsp_filter_indices, χch_filter_indices))
end
