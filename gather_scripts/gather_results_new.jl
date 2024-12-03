using Pkg
Pkg.activate("/scratch/projects/hhp00048/codes/LadderDGA.jl")
using LadderDGA, JLD2
using DataFrames, CSV

prefix = ARGS[1]
out_fn = ARGS[2]

df = DataFrame(U=Float64[], β=Float64[], Nk=Int[], Nf=Int[], μ_DMFT=Float64[], n=Float64[], EKin_DMFT=Float64[], EPot_DMFT=Float64[],
               k_χm_max=Tuple{Float64,Float64}[], k_χd_max=Tuple{Float64,Float64}[], χm_max=Float64[], χd_max=Float64[], χm_AF=Float64[], χm_FM=Float64[], 
               k_χm_min=Tuple{Float64,Float64}[], k_χd_min=Tuple{Float64,Float64}[], χm_min=Float64[], χd_min=Float64[], 
               λm_m=Float64[], λm_d=Float64[], μ_m=Float64[], EKin_m=Float64[], EPot_p1_m=Float64[], EPot_p2_m=Float64[], PP_p1_m=Float64[], PP_p2_m=Float64[], 
               ki_Σ_m_0_N=Int[], k_Σ_m_0_N=Tuple{Float64,Float64}[], ki_Σ_m_0_AN=Int[], k_Σ_m_0_AN=Tuple{Float64,Float64}[], Σ_m_0_N=ComplexF64[], Σ_m_0_AN=ComplexF64[], #converged_dm=Bool[], converged_sc_dm=Bool[],
               ki_Σ_m_lf_N=Int[], k_Σ_m_lf_N=Tuple{Float64,Float64}[], ki_Σ_m_lf_AN=Int[], k_Σ_m_lf_AN=Tuple{Float64,Float64}[], Σ_m_lf_N=ComplexF64[], Σ_m_lf_AN=ComplexF64[], conn_m=Float64[], #converged_dm=Bool[], converged_sc_dm=Bool[],
               λdm_m=Float64[], λdm_d=Float64[], μ_dm=Float64[], EKin_dm=Float64[], EPot_p1_dm=Float64[], EPot_p2_dm=Float64[], PP_p1_dm=Float64[], PP_p2_dm=Float64[], 
               ki_Σ_dm_0_N=Int[], k_Σ_dm_0_N=Tuple{Float64,Float64}[], ki_Σ_dm_0_AN=Int[], k_Σ_dm_0_AN=Tuple{Float64,Float64}[], Σ_dm_0_N=ComplexF64[], Σ_dm_0_AN=ComplexF64[], #converged_dm=Bool[], converged_sc_dm=Bool[],
               ki_Σ_dm_lf_N=Int[], k_Σ_dm_lf_N=Tuple{Float64,Float64}[], ki_Σ_dm_lf_AN=Int[], k_Σ_dm_lf_AN=Tuple{Float64,Float64}[], Σ_dm_lf_N=ComplexF64[], Σ_dm_lf_AN=ComplexF64[], conn_dm=Float64[], #converged_dm=Bool[], converged_sc_dm=Bool[],
               λdm_sc_m=Float64[], λdm_sc_d=Float64[], μ_dm_sc=Float64[], EKin_dm_sc=Float64[], EPot_p1_dm_sc=Float64[], EPot_p2_dm_sc=Float64[], PP_p1_dm_sc=Float64[], PP_p2_dm_sc=Float64[], 
               ki_Σ_dm_sc_0_N=Int[], k_Σ_dm_sc_0_N=Tuple{Float64,Float64}[], ki_Σ_dm_sc_0_AN=Int[], k_Σ_dm_sc_0_AN=Tuple{Float64,Float64}[], Σ_dm_sc_0_N=ComplexF64[], Σ_dm_sc_0_AN=ComplexF64[], #, #converged_dm_sc=Bool[], converged_sc_dm_sc=Bool[],
               ki_Σ_dm_sc_lf_N=Int[], k_Σ_dm_sc_lf_N=Tuple{Float64,Float64}[], ki_Σ_dm_sc_lf_AN=Int[], k_Σ_dm_sc_lf_AN=Tuple{Float64,Float64}[], Σ_dm_sc_lf_N=ComplexF64[], Σ_dm_sc_lf_AN=ComplexF64[], conn_dm_sc=Float64[]) #, #converged_dm_sc=Bool[], converged_sc_dm_sc=Bool[],
               # λdm_tsc_m=Float64[], λdm_tsc_d=Float64[], μ_dm_tsc=Float64[], EKin_dm_tsc=Float64[], EPot_p1_dm_tsc=Float64[], EPot_p2_dm_tsc=Float64[], PP_p1_dm_tsc=Float64[], PP_p2_dm_tsc=Float64[], Σ_dm_tsc_0_N=ComplexF64[], Σ_dm_tsc_0_AN=ComplexF64[] #converged_dm_tsc=Bool[], converged_sc_dm_tsc=Bool[]
               #  )

function find_k_N_AN(Σ_ladder, efmask, kG, mP) #LadderDGA.lin_fit
    tmp    =  Σ_ladder[:,0] .- (Inf + Inf*im) .* (.! efmask)
    Σ_AN, ki_AN  = findmin(abs.(imag(tmp)))
    tmp[:]  =  Σ_ladder[:,0] .+ (Inf + Inf*im) .* (.! efmask)
    Σ_N, ki_N  = findmin(imag(tmp))
    return Σ_ladder[ki_N,0], ki_N, Σ_ladder[ki_AN,0], ki_AN
end


# "res_ldga.jld2"
for (root, dirs, files) in walkdir(".")
    if endswith(root, "lDGA_julia")
        filelist = filter(f-> startswith(f, prefix) && endswith(f, ".jld2"), files)
        for ldga_fn in filelist
            println("adding dir: ", root)
            ki_N = 0
            ki_AN = 0
            row = jldopen(joinpath(root, ldga_fn), "r") do f
                h   = f["lDGAHelper"]
                nh = ceil(Int,size(f["χd"],2)/2)
                Nk = h.kG.Ns
                Nf = h.sP.n_iω
                χd_max, ki_d_max = findmax(f["χd"][:,nh])
                χm_max, ki_m_max = findmax(f["χm"][:,nh])
                χd_inv_min, ki_d_min = findmin(1 ./ f["χd"][:,nh])
                χm_inv_min, ki_m_min = findmin(1 ./ f["χm"][:,nh])

                k_χd_max = h.kG.kGrid[ki_d_max]
                k_χm_max = h.kG.kGrid[ki_m_max]
                k_χd_min = h.kG.kGrid[ki_d_min]
                k_χm_min = h.kG.kGrid[ki_m_min]

                res_m = f["res_m"]
                res_dm = f["res_dm"]
                res_dm_sc = f["res_dm_sc"]
                mode=LadderDGA.zero_freq
                Σ_m_N, ki_m_N, Σ_m_AN, ki_m_AN, Σ_m_N_lf, ki_m_N_lf, Σ_m_AN_lf, ki_m_AN_lf, rc_m = if !isnothing(res_m.Σ_ladder) 
                    efmask, rc = LadderDGA.estimate_connected_ef(res_m.Σ_ladder, h.kG, h.mP; ν0_estimator=mode)
                    Σ_m_N, ki_m_N, Σ_m_AN, ki_m_AN = find_k_N_AN(res_m.Σ_ladder, efmask, h.kG, h.mP)
                    Σ_m_N_lf, ki_m_N_lf, Σ_m_AN_lf, ki_m_AN_lf = find_k_N_AN(res_m.Σ_ladder, efmask, h.kG, h.mP)
                    Σ_m_N, ki_m_N, Σ_m_AN, ki_m_AN, Σ_m_N_lf, ki_m_N_lf, Σ_m_AN_lf, ki_m_AN_lf, rc
                else
                    (NaN, 0, NaN, 0, NaN, 0, NaN, 0, NaN)
                end
                Σ_dm_N, ki_dm_N, Σ_dm_AN, ki_dm_AN, Σ_dm_N_lf, ki_dm_N_lf, Σ_dm_AN_lf, ki_dm_AN_lf, rc_dm = if !isnothing(res_dm.Σ_ladder) 
                    efmask, rc = LadderDGA.estimate_connected_ef(res_dm.Σ_ladder, h.kG, h.mP; ν0_estimator=mode)
                    Σ_N, ki_N, Σ_AN, ki_AN = find_k_N_AN(res_dm.Σ_ladder, efmask, h.kG, h.mP)
                    Σ_N_lf, ki_N_lf, Σ_AN_lf, ki_AN_lf = find_k_N_AN(res_dm.Σ_ladder, efmask, h.kG, h.mP)
                    Σ_N, ki_N, Σ_AN, ki_AN, Σ_N_lf, ki_N_lf, Σ_AN_lf, ki_AN_lf, rc
                else
                    (NaN, 0, NaN, 0, NaN, 0, NaN, 0, NaN)
                end
                Σ_dmsc_N, ki_dmsc_N, Σ_dmsc_AN, ki_dmsc_AN, Σ_dmsc_N_lf, ki_dmsc_N_lf, Σ_dmsc_AN_lf, ki_dmsc_AN_lf, rc_dmsc = if !isnothing(res_dm_sc.Σ_ladder) 
                    efmask, rc = LadderDGA.estimate_connected_ef(res_dm_sc.Σ_ladder, h.kG, h.mP; ν0_estimator=mode)
                    Σ_N, ki_N, Σ_AN, ki_AN = find_k_N_AN(res_dm_sc.Σ_ladder, efmask, h.kG, h.mP)
                    Σ_N_lf, ki_N_lf, Σ_AN_lf, ki_AN_lf = find_k_N_AN(res_dm_sc.Σ_ladder, efmask, h.kG, h.mP)
                    Σ_N, ki_N, Σ_AN, ki_AN, Σ_N_lf, ki_N_lf, Σ_AN_lf, ki_AN_lf, rc
                else
                    (NaN, 0, NaN, 0, NaN, 0, NaN, 0, NaN)
                end

                println(Σ_m_N .- Σ_m_N_lf)
                k_Σ_m_0_N = h.kG.kGrid[ki_m_N]
                k_Σ_m_0_AN = h.kG.kGrid[ki_m_AN]
                k_Σ_dm_0_N = h.kG.kGrid[ki_dm_N]
                k_Σ_dm_0_AN = h.kG.kGrid[ki_dm_AN]
                k_Σ_dmsc_0_N = h.kG.kGrid[ki_dmsc_N]
                k_Σ_dmsc_0_AN = h.kG.kGrid[ki_dmsc_AN]
                k_Σ_m_N_lf = h.kG.kGrid[ki_m_N_lf]
                k_Σ_m_AN_lf = h.kG.kGrid[ki_m_AN_lf]
                k_Σ_dm_N_lf = h.kG.kGrid[ki_dm_N_lf]
                k_Σ_dm_AN_lf = h.kG.kGrid[ki_dm_AN_lf]
                k_Σ_dmsc_N_lf = h.kG.kGrid[ki_dmsc_N_lf]
                k_Σ_dmsc_AN_lf = h.kG.kGrid[ki_dmsc_AN_lf]

                println("nh = $nh, max_d = $(ki_d_max), max_m = $(ki_m_max) // $(k_χd_max) : $(k_χm_max)")
                χm_AF    = f["χm"][end,nh]
                χm_FM    = f["χm"][1,nh]
                row_ldga = [h.mP.U, h.mP.β, Nk, Nf, h.mP.μ, h.mP.n, h.mP.Ekin_DMFT, h.mP.Epot_DMFT, k_χm_max, k_χd_max, χm_max, χd_max, 
                            χm_AF, χm_FM, k_χm_min, k_χd_min, 1 / χm_inv_min, 1 / χd_inv_min]
                row_m = isnothing(res_m) ? repeat([NaN], 13) : [res_m.λm, res_m.λd, res_m.μ, res_m.EKin, res_m.EPot_p1, res_m.EPot_p2, res_m.PP_p1, res_m.PP_p2, 
                                                                ki_m_N, k_Σ_m_0_N, ki_m_AN, k_Σ_m_0_AN,  Σ_m_N, Σ_m_AN, ki_m_N_lf, k_Σ_m_N_lf, ki_m_AN_lf, k_Σ_m_AN_lf,  Σ_m_N_lf, Σ_m_AN_lf, rc_m]
                row_dm = isnothing(res_dm) ? repeat([NaN], 13) : [res_dm.λm, res_dm.λd, res_dm.μ, res_dm.EKin, res_dm.EPot_p1, res_dm.EPot_p2, res_dm.PP_p1, res_dm.PP_p2, 
                                                                  ki_dm_N,  k_Σ_dm_0_N, ki_dm_AN,  k_Σ_dm_0_AN, Σ_dm_N, Σ_dm_AN, ki_dm_N_lf,  k_Σ_dm_N_lf, ki_dm_AN_lf,  k_Σ_dm_AN_lf, Σ_dm_N_lf, Σ_dm_AN_lf, rc_dm]
                row_dm_sc = isnothing(res_dm_sc) ? repeat([NaN], 13) : [res_dm_sc.λm, res_dm_sc.λd, res_dm_sc.μ, res_dm_sc.EKin, res_dm_sc.EPot_p1, res_dm_sc.EPot_p2, res_dm_sc.PP_p1, res_dm_sc.PP_p2, 
                                                                        ki_dmsc_N,  k_Σ_dmsc_0_N, ki_dmsc_AN,  k_Σ_dmsc_0_AN, Σ_dmsc_N, Σ_dmsc_AN, ki_dmsc_N_lf,  k_Σ_dmsc_N_lf, ki_dmsc_AN_lf,  k_Σ_dmsc_AN_lf, Σ_dmsc_N, Σ_dmsc_AN_lf, rc_dmsc]
                vcat(row_ldga,row_m,row_dm,row_dm_sc)
            end
            push!(df, row)
        end
    end
end
println(df)
CSV.write(out_fn, df)

