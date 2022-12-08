using Pkg
Pkg.activate("/scratch/projects/hhp00048/codes/lDGAPostprocessing")
using lDGAPostprocessing

using DelimitedFiles
using LinearAlgebra

input_str = ARGS[1]
U = tryparse(Float64,ARGS[2])

ϵₖ, Vₖ, μ    = read_anderson_parameters("./hubb.andpar");
#nh = floor(Int, length(iν_in)/2 + 1)

function model_ED(iν, p)
    N_andpar = Int(length(p)/2)
    Δ_fit = Array{ComplexF64,length(iν)}
    GW_fit = Array{ComplexF64,length(iν)}
    for (i,νn) in enumerate(iν)
        Δ_fit[i]  = sum(p[N_andpar+1:end] .* conj(p[N_andpar+1:end]) ./ (νn .- p[1:N_andpar]))
        GW_fit[i] = 1 / (νn + μ - Δ_fit[i])
    end
    return GW_fit, Δ_fit
end
