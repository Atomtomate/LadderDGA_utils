using LsqFit
using Optim
using DelimitedFiles
using LinearAlgebra
using HDF5

if length(ARGS) == 0
    println("Call with arguments [file name] [β] [U] [μ]")
    return
end

input_str = ARGS[1]

β = tryparse(Float64,ARGS[2])
U = tryparse(Float64,ARGS[3])
μ = tryparse(Float64,ARGS[4])
#n = tryparse(Float64,ARGS[5])


function calc_E_ED(iνₙ, ϵₖ, Vₖ, GImp, U, n, μ, β)
    E_kin = 0.0
    E_pot = 0.0
    vk = sum(Vₖ .^ 2)
    Σ_hartree = n * U/2
    E_pot_tail = (U^2)/2 * n * (1-n/2) - Σ_hartree*(Σ_hartree-μ)
    E_kin_tail = vk

    for n in 1:length(GImp)
        Δ_n = sum((Vₖ .^ 2) ./ (iνₙ[n] .- ϵₖ))
        Σ_n = iνₙ[n] .- Δ_n .- 1.0 ./ GImp[n] .+ μ
        E_kin += 2*real(GImp[n] * Δ_n - E_kin_tail/(iνₙ[n]^2))
        E_pot += 2*real(GImp[n] * Σ_n - E_pot_tail/(iνₙ[n]^2))
    end
    E_kin = E_kin .* (2/β) - (β/2) .* E_kin_tail
    E_pot = E_pot .* (1/β) .+ 0.5*Σ_hartree .- (β/4) .* E_pot_tail
    return E_kin, E_pot
end

iν_in, g, g_imp = if endswith(input_str, ".dat")
    g_in     = readdlm(g_in_f)
    iν_in, g, g_imp = 1im * g_in[:,1], g_in[:,2] .+ 1im * g_in[:,3], nothing;
elseif endswith(input_str, ".hdf5") 
    f    = h5open(ARGS[1])
    g_in = read(f["dmft-last/ineq-001/giw/value"])
    g_imp_in = read(f["dmft-last/ineq-001/giw/value"])
    g    = (g_in[:,1,1] .+ g_in[:,2,1] ) ./ 2
    g_imp    = (g_imp_in[:,1,1] .+ g_imp_in[:,2,1] ) ./ 2
    iν_in = read(f[".axes/iw"])
    1im * iν_in, g, g_imp
else
    error("File extension not recognized")
end



nh = floor(Int, length(iν_in)/2 + 1)
Δ        = -(1 ./ g .- iν_in .- μ);
model_ED_hf(iν, p) = - 2 .* iν .* sum(p[3:4] .^ 2 ./ (transpose(iν .^ 2) .+ p[1:2] .^ 2), dims=1)[1,:]
best_N = 0
best_check = Inf
best_ϵp = nothing
best_Vp = nothing
best_range = nothing

for Nν_i in 1:100
    range = (nh-Nν_i):(nh+Nν_i-1)
    Δ_i  = Δ[range]
    ν_i = iν_in[range]

    fit3 = curve_fit(model_ED_hf, imag.(ν_i), imag.(Δ_i), [0.5, 1.0, 0.25, 0.25])
    ϵp = round.([fit3.param[1], fit3.param[2], -fit3.param[1], -fit3.param[2]], digits=15)
    Vp = round.([fit3.param[3], fit3.param[4], fit3.param[3], fit3.param[4]], digits=15)
    check = abs(sum(Vp .^ 2) - 0.25)
    if check < best_check
        println("=> New best")
        global best_N = Nν_i
        global best_check = check
        global best_ϵp = ϵp
        global best_Vp = Vp
        global best_range = range 
    end
    println("N = $Nν_i: sum(V^2_k) - 0.25 = $(round(check,digits=6))")
end
println("    Anderson Parameter Checks: ")
println("   ============================   ")
println("1. min(|V_k|)       = $(round(minimum(abs.(best_Vp)), digits=4))")
println("2. sum(V^2_k)       = $(round(sum(best_Vp .^ 2), digits=4))")
println("3. min(|e_k|)       = $(round(minimum(abs.(best_ϵp)), digits=4))")
println("4. min(|e_i - e_j|) = $(round(minimum(abs.(best_ϵp .- best_ϵp')  + Inf .* I ), digits=4))")
println("   ============================   ")


out = """
           ========================================
               1-band            30-Sep-95 LANCZOS
            ========================================
NSITE     5 IWMAX32768
 $(β)d0, -12.0, 12.0, 0.007
c ns,imaxmu,deltamu, # iterations, conv.param.
  5, 0, 0.d0, 80,  1.d-14
c ifix(0,1), <n>,   inew, iauto
Eps(k)
   $(best_ϵp[1])
   $(best_ϵp[2])
   $(best_ϵp[3])
   $(best_ϵp[4])
 tpar(k)
   $(best_Vp[1])
   $(best_Vp[2])
   $(best_Vp[3])
   $(best_Vp[4])
  $μ          #chemical potential
"""
open("hubb.andpar","w") do io 
    print(io,out)
end

#EKin, EPot = calc_E_ED(iν_in[best_range], best_ϵp, best_Vp, g_imp[best_range], U, n, μ, β)

println("Anderson Parameters: \n", best_ϵp, "\n", best_Vp)
#println("Kinetic Energy: $(EKin)\nPotential Energy: $(EPot)")
# run_dir = pwd()
# include("../ed_consistency.jl")
# println("Anderson Parameters: \n")
# println("ϵp: \n", best_ϵp)
# println("Vp: \n", best_Vp)
