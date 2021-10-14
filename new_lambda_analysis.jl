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

function find_epot(λch_vals, c2_curve, res)
    rrr = []
    yvals = c2_curve
    xvals = λch_vals
    sc_ind = findfirst(x->x<0,sign.(yvals))
    (sc_ind == nothing) && return [NaN, NaN, NaN, NaN]
    (sc_ind == 1) && return [NaN, NaN, NaN, NaN]
    return res[sc_ind-1,5], res[sc_ind-1,6], res[sc_ind,5], res[sc_ind,6]
end

function λsp_of_λch(nlQ_sp::NonLocalQuantities, nlQ_ch::NonLocalQuantities, kG, mP, sP; max_λsp=30.0, λch_max=1.0, n_λch=20)
    ωindices = (sP.dbg_full_eom_omega) ? (1:size(nlQ_sp.χ,2)) : intersect(nlQ_sp.usable_ω, nlQ_ch.usable_ω)
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
                        Gνω::GνqT, FUpDo::AbstractArray{Complex{Float64},3}, kG::ReducedKGrid,
                        mP::ModelParameters, sP::SimulationParameters)
    println("starting $(mP)")
    println(stderr,"starting $(mP)")
    flush(stdout)
    flush(stderr)
    # Prepare data
    res = Array{Float64, 2}(undef, length(λch_range), 6)
    ωindices = (sP.dbg_full_eom_omega) ? (1:size(bubble,3)) : intersect(nlQ_sp.usable_ω, nlQ_ch.usable_ω)
    lur = length(ωindices)
    νmax = floor(Int,lur/3)
    if νmax > 5
        # --- prepare auxiliary vars ---
        Kνωq = Array{ComplexF64, length(gridshape(kG))}(undef, gridshape(kG)...)
        Kνωq_pre = Array{ComplexF64, 1}(undef, size(bubble,LadderDGA.q_axis))
        νGrid = 0:(νmax-1)
        iν_n = LadderDGA.iν_array(mP.β, νGrid)
        iωn = 1im .* 2 .* (-sP.n_iω:sP.n_iω)[ωindices] .* π ./ mP.β
        iωn2_sub = real.([i == 0 ? 0 : mP.Ekin_DMFT / (i)^2 for i in iωn])
        nd = length(gridshape(kG))
        
        Σ_ladder_i = Array{Complex{Float64},2}(undef, size(bubble,1), νmax)
        χsp_λ = similar(nlQ_sp.χ[:,ωindices])
        χch_λ = similar(nlQ_ch.χ[:,ωindices])
    
        # Prepare data
        corr = Σ_correction(ωindices, bubble, FUpDo, sP)
        nh    = ceil(Int64, size(nlQ_sp.χ,2)/2)
        χsp_min    = -1 / maximum(real.(nlQ_sp.χ[:,nh]))
        χch_min    = -1 / maximum(real.(nlQ_ch.χ[:,nh]))

        Σ_hartree = mP.n * mP.U/2
        E_pot_tail_c = [zeros(size(kG.ϵkGrid)),
                        (mP.U^2 * 0.5 * mP.n * (1-0.5*mP.n) .+ Σ_hartree .* (kG.ϵkGrid .+ Σ_hartree .- mP.μ))]
        tail = [1 ./ (iν_n .^ n) for n in 1:length(E_pot_tail_c)]
        E_pot_tail = permutedims(sum((E_pot_tail_c[i])' .* tail[i] for i in 1:length(tail)),(2,1))
        E_pot_tail_inv = sum((mP.β/2)  .* [Σ_hartree .* ones(size(kG.ϵkGrid)), (-mP.β/2) .* E_pot_tail_c[2]])

        Σ_corr = Σ_loc[1:length(Σ_ladderLoc)] .- Σ_ladderLoc[:] .+ Σ_hartree

        for i in 1:length(λch_range)
            λ = [spOfch[i], λch_range[i]]
            fill!(Σ_ladder_i, zero(eltype(Σ_ladder_i)))
            χ_λ!(χsp_λ, view(nlQ_sp.χ,:,ωindices), λ[1])
            χ_λ!(χch_λ, view(nlQ_ch.χ,:,ωindices), λ[2])
            lhs_c1, lhs_c2 = 0.0, 0.0
            for ωii in 1:length(ωindices)
                ωi = ωindices[ωii]
                ωn = (ωi - sP.n_iω) - 1
                fsp = 1.5 .* (1 .+ mP.U .* view(χsp_λ,:,ωii))
                fch = 0.5 .* (1 .- mP.U .* view(χch_λ,:,ωii))
                νZero = LadderDGA.ν0Index_of_ωIndex(ωi, sP)
                maxn = minimum([νZero + νmax - 1, size(nlQ_ch.γ, LadderDGA.ν_axis)])
                for (νi,νn) in enumerate(νZero:maxn)
                    v = selectdim(Gνω,nd+1,(νi-1) + ωn + sP.fft_offset)
                    @simd for qi in 1:length(Kνωq_pre)
                        @inbounds Kνωq_pre[qi] = nlQ_sp.γ[qi,νn,ωi] * fsp[qi] - nlQ_ch.γ[qi,νn,ωi] * fch[qi] - 1.5 + 0.5 + corr[qi,νn,ωii]
                    end
                    expandKArr!(kG,Kνωq,Kνωq_pre)
                    Dispersions.mul!(Kνωq, kG.fftw_plan, Kνωq)
                    @simd for ki in 1:length(Kνωq)
                        @inbounds Kνωq[ki] *= v[ki]
                    end
                    Dispersions.ldiv!(Kνωq, kG.fftw_plan, Kνωq)
                    reduceKArr!(kG, Kνωq_pre, Dispersions.ifft_post(kG, Kνωq)) 
                    @simd for i in 1:length(Kνωq_pre)
                        @inbounds Σ_ladder_i[i,νi] += mP.U * Kνωq_pre[i]/ (kG.Nk*mP.β)
                    end
                end
                tsp, tch = 0.0, 0.0
                for qi in 1:length(kG.kMult)
                    tsp += kG.kMult[qi]*real(χsp_λ[qi,ωii])
                    tch += kG.kMult[qi]*real(χch_λ[qi,ωii])
                end
                lhs_c1 += (tch + tsp) / (2*kG.Nk) - iωn2_sub[ωii]
                lhs_c2 += (tch - tsp) / (2*kG.Nk)
            end
            E_pot = 0.0
            for qi in 1:length(kG.kMult)
                GΣ_λ = 0.0
                for i in 1:νmax
                    Σ_ladder_i[qi,i] += Σ_corr[i]
                    GΣ_λ += 2 * real(Σ_ladder_i[qi,i] * LadderDGA.G_from_Σ(iν_n[i], mP.β, mP.μ, kG.ϵkGrid[qi], Σ_ladder_i[qi,i]) - E_pot_tail[qi,i])

                end
                GΣ_λ += E_pot_tail_inv[qi]   # ν summation
                E_pot += kG.kMult[qi]*GΣ_λ   # k intgration
            end
            E_pot = E_pot / (kG.Nk * mP.β)
            lhs_c1 = lhs_c1/mP.β - mP.Ekin_DMFT*mP.β/12
            lhs_c2 = lhs_c2/mP.β
            rhs_c1 = mP.n/2 * (1 - mP.n/2)
            rhs_c2 = E_pot/mP.U - (mP.n/2) * (mP.n/2)
            res[i,:] = [λ[1] λ[2] lhs_c1 rhs_c1 lhs_c2 rhs_c2]
        end
    end
    return res
end

function new_λ_from_c2(c2_res, impQ_sp, impQ_ch, FUpDo, Σ_loc, Σ_ladder_loc, nlQ_sp, nlQ_ch_λ, bubble, gLoc_fft, kG, mP, sP)
    λsp, λch = if size(c2_res,1) >= 1
        λch = find_zero(c2_res[:,2], c2_res[:,5] .- c2_res[:,6])
        nlQ_ch_λ = deepcopy(nlQ_ch)
        nlQ_ch_λ.χ = LadderDGA.χ_λ(nlQ_ch.χ, λch)
        nlQ_ch_λ.λ = λch
        λsp = LadderDGA.λ_correction(:sp,impQ_sp, impQ_ch, FUpDo, Σ_loc, Σ_ladder_loc, nlQ_sp, nlQ_ch_λ, bubble, gLoc_fft, kG, mP, sP)
    else
        NaN, NaN
    end
    λsp, λch
end
