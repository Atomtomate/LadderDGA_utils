using LsqFit
using DelimitedFiles
using LinearAlgebra
using HDF5

input_str = ARGS[1]

β = tryparse(Float64,ARGS[2])
U = tryparse(Float64,ARGS[3])
Nν = tryparse(Int,ARGS[4])
μ = U/2

iν_in, g = if endswith(input_str, ".dat")
    g_in     = readdlm(g_in_f)
    iν_in, g = 1im * g_in[:,1], g_in[:,2] .+ 1im * g_in[:,3];
elseif endswith(input_str, ".hdf5") 
    f    = h5open(ARGS[1])
    g_in = read(f["dmft-last/ineq-001/giw/value"])
    g    = (g_in[:,1,1] .+ g_in[:,2,1] ) ./ 2
    iν_in = read(f[".axes/iw"])
    1im * iν_in, g
else
    error("File extension not recognized")
end
nh = floor(Int, length(iν_in)/2 + 1)
Δ        = -(1 ./ g .- iν_in .- μ);
model_ED3(iν, p) = - 2 .* iν .* sum(p[3:4] .^ 2 ./ (transpose(iν .^ 2) .+ p[1:2] .^ 2), dims=1)[1,:]
best_N = 0
best_check = Inf
best_ϵp = nothing
best_Vp = nothing

for Nν_i in 1:nh-2
    range = (nh-Nν_i):(nh+Nν_i-1)
    Δ_i  = Δ[range]
    ν_i = iν_in[range]

    fit3 = curve_fit(model_ED3, imag.(ν_i), imag.(Δ_i), [0.25, 0.25, -1.3, -0.5,])
    ϵp = round.([fit3.param[1], fit3.param[2], -fit3.param[1], -fit3.param[2]], digits=15)
    Vp = round.([fit3.param[3], fit3.param[4], fit3.param[3], fit3.param[4]], digits=15)
    check = abs(sum(Vp .^ 2) - 0.25)
    if check < best_check
        println("=> New best")
        global best_N = Nν_i
        global best_check = check
        global best_ϵp = ϵp
        global best_Vp = Vp
    end
    println("N = $Nν_i: sum(V^2_k) - 0.25 = $(round(check,digits=4))")
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
  $β d0, -12.0, 12.0, 0.007
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


println("Anderson Parameters: \n")
println("ϵp: \n", best_ϵp)
println("Vp: \n", best_Vp)
