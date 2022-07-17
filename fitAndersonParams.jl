using LsqFit
using DelimitedFiles
using HDF5

input_str = ARGS[1]

β = tryparse(Float64,ARGS[2])
U = tryparse(Float64,ARGS[3])
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
Δ        = -(1 ./ g .- iν_in .- μ);
Δ, iν_in

model_ED3(iν, p) = - 2 .* iν .* sum(p[3:4] .^ 2 ./ (transpose(iν .^ 2) .+ p[1:2] .^ 2), dims=1)[1,:]
fit3 = curve_fit(model_ED3, imag.(iν_in), imag.(Δ), [0.25, 0.25, -1.3, -0.5,])
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


println("Anderson Parameters: \n", fit3.param)
