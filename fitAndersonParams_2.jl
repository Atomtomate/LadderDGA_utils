using LsqFit
# using Optim
using DelimitedFiles
using LinearAlgebra
using HDF5

input_str = ARGS[1]
β = tryparse(Float64,ARGS[2])
U = tryparse(Float64,ARGS[3])
μ = tryparse(Float64,ARGS[4])
N = tryparse(Int, ARGS[5])
outf = ARGS[6]


iν_in, g, g_imp, Δ_in  = if endswith(input_str, ".dat")
    g_in     = readdlm(g_in_f)
    iν_in, g, g_imp = 1im * g_in[:,1], g_in[:,2] .+ 1im * g_in[:,3], nothing;
elseif endswith(input_str, ".hdf5") 
    f    = h5open(ARGS[1])
    g_in = read(f["dmft-last/ineq-001/g0iw/value"])
    g_imp_in = read(f["dmft-last/ineq-001/g0iw/value"])
    #g    = (g_in[:,1,1]) 
    #g_imp    = (g_imp_in[:,1,1]) 
    Δ_in_tmp = read(f["dmft-last/ineq-001/fiw/value"])
    Δ_in = (Δ_in_tmp[:,1,1] .+ Δ_in_tmp[:,2,1]) ./ 2
    g    = (g_in[:,1,1] .+ g_in[:,2,1] ) ./ 2
    g_imp    = (g_imp_in[:,1,1] .+ g_imp_in[:,2,1] ) ./ 2
    iν_in = read(f[".axes/iw"])
    1im * iν_in, g, g_imp, Δ_in
else
    error("File extension not recognized")
end

nh = floor(Int, length(iν_in)/2 + 1)
Δ        =  conj.(-1 ./ g .+ iν_in .+ μ);
#println(sum(abs.(Δ[950:1050] .- Δ_in[950:1050])))
#Δ = Δ_in
function model_ED(iν, p)
    Δ_fit = zeros(ComplexF64,length(iν))
    for (i,νn) in enumerate(iν)
        tmp  = sum((p[(N+1):end] .^ 2) ./ (νn .- p[1:N]))
        Δ_fit[i] = tmp
    end
    return conj.(Δ_fit)
end

function model_ED_real(iν, p)
    Δ_fit = zeros(ComplexF64,length(iν))
    for (i,νn) in enumerate(iν)
        tmp  = sum((p[(N+1):end] .^ 2) ./ (νn .- p[1:N]))
        Δ_fit[i] = tmp
    end
    return vcat(real(Δ_fit), -imag(Δ_fit))
end
#model_ED_hf(iν, p) = - 2 .* iν .* sum(p[3:4] .^ 2 ./ (transpose(iν .^ 2) .+ p[1:2] .^ 2), dims=1)[1,:]
best_N = 0
best_check = Inf
best_fit = nothing
best_ϵp = nothing
best_Vp = nothing
best_range = nothing


# for Nν_i in 100:300
Nν_i = 100
    νrange = (nh-Nν_i):(nh+Nν_i-1)
    Δ_i  = Δ[νrange]
    ν_i = iν_in[νrange]
    p0 = vcat(range(-U/2,U/2,length=N),range(0.1,1.0,length=N))
    #fit3 = curve_fit(model_ED_hf, imag.(ν_i), imag.(Δ_i), [0.25, 0.25, -1.3, -0.5,])
    fit = curve_fit(model_ED_real, ν_i, vcat(real(Δ_i),imag(Δ_i)), p0)
    ϵp = round.(fit.param[1:N], digits=15)
    Vp = round.(fit.param[(N+1):end], digits=15)
    # ϵp = round.([fit3.param[1], fit3.param[2], -fit3.param[1], -fit3.param[2]], digits=15)
    # Vp = round.([fit3.param[3], fit3.param[4], fit3.param[3], fit3.param[4]], digits=15)
    check = abs(sum(Vp .^ 2) - 0.25)
    if check < best_check
        println("=> New best")
        global best_N = Nν_i
        global best_check = check
        global best_ϵp = ϵp
        global best_Vp = Vp
        global best_range = νrange 
        global best_fit = fit
    end
    println("N = $Nν_i: sum(V^2_k) - 0.25 = $(round(check,digits=6))")
# end

println("    Anderson Parameter Checks: ")
println("   ============================   ")
println("1. min(|V_k|)       = $(round(minimum(abs.(best_Vp)), digits=4))")
println("2. sum(V^2_k)       = $(round(sum(best_Vp .^ 2), digits=4))")
println("3. min(|e_k|)       = $(round(minimum(abs.(best_ϵp)), digits=4))")
println("4. min(|e_i - e_j|) = $(round(minimum(abs.(best_ϵp .- best_ϵp')  + Inf .* I ), digits=4))")
println("   ============================   ")

function write_andpar(outf, best_ϵp, best_Vp)
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
    """
    for ϵki in best_ϵp
        out = out * "$ϵki \n"
    end
    out = out * "tpar(k)\n"
    for Vki in best_Vp
        out = out * "$Vki \n"
    end
    out = out * "$μ"
    open(outf,"w") do io 
        print(io,out)
    end
end


println("Anderson Parameters: \n", best_ϵp, "\n", best_Vp)
write_andpar(outf, best_ϵp, best_Vp)
# run_dir = pwd()
# include("../ed_consistency.jl")
# println("Anderson Parameters: \n")
# println("ϵp: \n", best_ϵp)
# println("Vp: \n", best_Vp)
