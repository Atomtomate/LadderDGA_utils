using Dispersions
using LadderDGA

include((@__DIR__)*"/helpers/new_lambda_helpers.jl")

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

function λsp_of_λch(nlQ_sp::NonLocalQuantities, nlQ_ch::NonLocalQuantities, kG, mP, sP; max_λsp=30.0, λch_max=1.0, n_λch=20)
    ωindices = (sP.dbg_full_eom_omega) ? (1:size(bubble,3)) : intersect(nlQ_sp.usable_ω, nlQ_ch.usable_ω)
    iωn = 1im .* 2 .* (-sP.n_iω:sP.n_iω)[ωindices] .* π ./ mP.β
    nh  = ceil(Int64, size(nlQ_sp.χ,2)/2)
    χsp_min    = -minimum(1 ./ real.(nlQ_sp.χ[:,nh]))
    χch_min    = -minimum(1 ./ real.(nlQ_ch.χ[:,nh]))

    λch_range = sort(union(range(χch_min,λch_max,n_λch), [0]))
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
                        nlQ_sp::NonLocalQuantities, nlQ_ch::NonLocalQuantities, bubble::BubbleT,
                        Σ_ladderLoc, Σ_loc,
                        GLoc_fft::GνqT, FUpDo::AbstractArray{Complex{Float64},3}, kG::ReducedKGrid,
                        mP::ModelParameters, sP::SimulationParameters)
    println("starting $(mP)")
    println(stderr,"starting $(mP)")
    flush(stdout)
    flush(stderr)
    # Prepare data
    res = Array{Float64, 2}(undef, length(λch_range), 6)
    ωindices = (sP.dbg_full_eom_omega) ? (1:size(bubble,3)) : intersect(nlQ_sp.usable_ω, nlQ_ch.usable_ω)
    iωn = 1im .* 2 .* (-sP.n_iω:sP.n_iω)[ωindices] .* π ./ mP.β
    corr = Σ_correction(ωindices, bubble, FUpDo, sP)
    (sP.tc_type_f != :nothing) && LadderDGA.extend_corr!(corr)

    νmax = trunc(Int,size(bubble,LadderDGA.ν_axis)/3)
    νGrid = 0:(νmax-1)
    iν_n = LadderDGA.iν_array(mP.β, νGrid)
    Σ_hartree = mP.n * mP.U/2

    E_pot_tail_c = [zeros(size(kG.ϵkGrid)),
                (mP.U^2 * 0.5 * mP.n * (1-0.5*mP.n) .+ Σ_hartree .* (kG.ϵkGrid .+ Σ_hartree .- mP.μ))]
    tail = [1 ./ (iν_n .^ n) for n in 1:length(E_pot_tail_c)]
    E_pot_tail =  sum(E_pot_tail_c[i] .* transpose(tail[i]) for i in 1:length(tail))
	E_pot_tail_inv = sum((mP.β/2)  .* [Σ_hartree .* ones(size(kG.ϵkGrid)), (-mP.β/2) .* E_pot_tail_c[2]])

    χsp = similar(nlQ_sp.χ);
    χch = similar(nlQ_ch.χ);
    χupdo_ω = Array{eltype(χsp),1}(undef, length(ωindices))
    χupup_ω = Array{eltype(χsp),1}(undef, length(ωindices))

    Σ_list = Array{Complex{Float64},3}(undef, length(λch_range), length(Σ_ladderLoc), length(kG.kGrid))
    G_λ = Array{ComplexF64, 2}(undef, length(kG.kMult), νmax)

    Σ_λ_ω = Array{Complex{Float64},3}(undef,size(bubble,1),νmax,size(bubble,3))
    Σ_λ = Array{Complex{Float64},2}(undef,size(bubble,1),νmax)
    Kνωq = Array{ComplexF64, length(gridshape(kG))}(undef, gridshape(kG)...)
    Kνωq_pre = Array{ComplexF64, 1}(undef, length(kG.kMult))

    for i in 1:length(λch_range)
        λsp_i = spOfch[i]
        λch_i = λch_range[i]
        LadderDGA.χ_λ!(χsp, nlQ_sp.χ, λsp_i, ωindices)
        LadderDGA.χ_λ!(χch, nlQ_ch.χ, λch_i, ωindices)
        #TODO: sc here, this could be inclines to avoid realloc of Σ_ladder_ω

        LadderDGA.calc_Σ_ω!(Σ_λ_ω, Kνωq, Kνωq_pre, ωindices, χsp, nlQ_sp.γ, χch, nlQ_ch.γ, GLoc_fft, corr, mP.U, kG, sP)
        @inbounds @views Σ_λ[:,:] = (mP.U/mP.β) .* sum(Σ_λ_ω, dims=[3])[:,:,1] .+ Σ_hartree
        Σ_λ = Σ_loc_correction(Σ_λ, Σ_ladderLoc, Σ_loc)

        #@simd for qi in 1:size(Σ_λ, 1)
        #    @inbounds view(Σ_λ, qi, :) = view(Σ_λ, qi, :) .+ Σ_corr
        #end
        for (j,ek) in enumerate(kG.ϵkGrid)
            for ni in 0:νmax-1
                @inbounds G_λ[j,ni+1] = LadderDGA.G_from_Σ(ni, mP.β, mP.μ, ek, Σ_λ[j,ni+1])
            end
        end
        println("aaa ", Σ_λ[4,2:8])
        G_λ2 = transpose(LadderDGA.flatten_2D(LadderDGA.G_from_Σ(Σ_λ, kG.ϵkGrid, νGrid, mP)));
        #TODO: sc end here
        EPot = LadderDGA.calc_E_pot(kG, G_λ, Σ_λ, E_pot_tail, E_pot_tail_inv, mP.β)
        for (wi,w) in enumerate(ωindices)
            @inbounds χupup_ω[wi] = real(kintegrate(kG, view(χch,:,w) .+ view(χsp,:,w))[1]) / 2
            @inbounds χupdo_ω[wi] = real(kintegrate(kG, view(χch,:,w) .- view(χsp,:,w))[1]) / 2
        end
        #χupup_ω = subtract_tail(χupup_ω, mP.Ekin_DMFT, iωn)
        lhs_c1 = sum(χupup_ω)/mP.β #- mP.Ekin_DMFT*mP.β/12
        lhs_c2 = sum(χupdo_ω)/mP.β
        rhs_c1 = mP.n/2 * (1 - mP.n/2)
        rhs_c2 = EPot/mP.U - (mP.n/2) * (mP.n/2)
        res[i,:] = [λsp_i λch_i lhs_c1 rhs_c1 lhs_c2 rhs_c2]
    end
    println("finishing $(mP)")
    return res
end
