using LsqFit
using DelimitedFiles

g_in_f = ARGS[1]

U = 3.3
μ = U/2
g_in     = readdlm(g_in_f)
iν_in, g = 1im * g_in[:,1], g_in[:,2] .+ 1im * g_in[:,3];
Δ         = -(1 ./ g .- iν_in .- μ);

model_ED3(iν, p) = - 2 .* iν .* sum(p[3:4] .^ 2 ./ (transpose(iν .^ 2) .+ p[1:2] .^ 2), dims=1)[1,:]
fit3 = curve_fit(model_ED3, imag.(iν_in), imag.(Δ), [0.25, 0.25, -1.3, -0.5,])
println("Anderson Parameters: \n", fit3.param)
