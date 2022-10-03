using LsqFit
using DelimitedFiles
using HDF5

input_str = ARGS[1]

β = tryparse(Float64,ARGS[2])
U = tryparse(Float64,ARGS[3])
N = tryparse(Int,ARGS[4])
μ = U/2
n = 1.0

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
    g_in = read(f["dmft-last/ineq-001/g0iw/value"])
    g_imp_in = read(f["dmft-last/ineq-001/g0iw/value"])
    g    = (g_in[:,1,1] .+ g_in[:,2,1] ) ./ 2
    g_imp    = (g_imp_in[:,1,1] .+ g_imp_in[:,2,1] ) ./ 2
    iν_in = read(f[".axes/iw"])
    1im * iν_in, g, g_imp
else
    error("File extension not recognized")
end

Δ        = -(1 ./ g .- iν_in .- μ);
nh = length(Δ)/2
slice = Int.((nh+1-N):(nh+N))
iν_fit = iν_in[slice]
Δ_fit = Δ[slice]
g_imp_fit = g_imp[slice]

model_ED(n,p) = sum(( p[5:8] .^ 2 ) ./ (transpose(n) .- p[1:4]), dims = 1)[1, :]
model_ED3(iν, p) = - 2 .* iν .* sum(p[3:4] .^ 2 ./ (transpose(iν .^ 2) .+ p[1:2] .^ 2), dims=1)[1,:]

fit3 = curve_fit(model_ED3, imag.(iν_fit), imag.(Δ_fit), [0.25, 0.25, -1.3, -0.5,])
fit4 = curve_fit(model_ED, imag.(iν_fit), imag.(Δ_fit), [0.25, 0.25, 0.25, 0.25, -1.3, -0.5, 1.3, 0.5])

s = sprint(show, fit3; context=:compact => false)
ϵₖ = [fit3.param[1], fit3.param[2], -fit3.param[1], -fit3.param[2]]
Vₖ = [fit3.param[3], fit3.param[4],  fit3.param[3],  fit3.param[4]]

println("check: $(2 * sum(fit3.param[3:4].^2))")
println("converged: $(fit3.converged)")
out = """
           ========================================
               1-band            30-Sep-95 LANCZOS
            ========================================
NSITE     5 IWMAX32768
  $β d0, -12.0, 12.0, 0.007
c ns,imaxmu,deltamu, # iterations, conv.param.
  5, 0, 0.d0, 80,  1.d-14
c ifix(0,1), <n>,   inew, iauto
Eps(k)
   $(round(fit3.param[1],digits=15))
   $(round(fit3.param[2],digits=15))
   $(-round(fit3.param[1],digits=15))
   $(-round(fit3.param[2],digits=15))
 tpar(k)
   $(round(fit3.param[3],digits=15))
   $(round(fit3.param[4],digits=15))
   $(round(fit3.param[3],digits=15))
   $(round(fit3.param[4],digits=15))
  $μ          #chemical potential
"""
open("hubb.andpar","w") do io 
    print(io,out)
end

EKin, EPot = calc_E_ED(iν_fit, ϵₖ, Vₖ, g_imp_fit, U, n, μ, β)

println("Anderson Parameters: \n", ϵₖ, "\n", Vₖ)
println("Kinetic Energy: $(EKin)\nPotential Energy: $(EPot)")
run_dir = pwd()
include("../ed_consistency.jl")
