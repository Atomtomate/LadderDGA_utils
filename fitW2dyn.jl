using LsqFit
using LinearAlgebra
using HDF5

f = h5open(ARGS[1])
N = tryparse(Int,ARGS[2])
Nν_max = tryparse(Int,ARGS[3])
outf = ARGS[4]

β = read_attribute(f[".config"], "general.beta")
U = read_attribute(f[".config"], "atoms.1.udd")
μ = read(f["dmft-last/mu/value"])
Δ_in_tmp = read(f["dmft-last/ineq-001/fiw/value"])
Δ_in = (Δ_in_tmp[:,1,1] .+ Δ_in_tmp[:,2,1]) ./ 2
iν_in = 1im .* read(f[".axes/iw"])

nh = floor(Int, length(iν_in)/2 + 1)

function model_ED(iν, p)
    Δ_fit = zeros(ComplexF64,length(iν))
    for (i,νn) in enumerate(iν)
        Δ_fit[i] = sum((p[(N+1):end] .^ 2) ./ (νn .- p[1:N]))
    end
    return conj.(Δ_fit)
end

function model_ED_real(iν, p)
    Δ_fit = model_ED(iν,p)
    return vcat(real(Δ_fit),imag(Δ_fit))
end

best_N = 0
best_check = Inf
best_fit = nothing
best_ϵp = nothing
best_Vp = nothing
best_range = nothing
Vp_min = 0.010

println("Trying to determine best fit for NBath = $N, U = $U, β=$β, min(abs(Vp)) = $Vp_min in range of fermionic frequencies $(10:(10+Nν_max))")

for Nν_i in 10:(10+Nν_max)
    νrange = (nh-Nν_i):(nh+Nν_i-1)
    Δ_i  = Δ_in[νrange]
    ν_i = iν_in[νrange]
    p0 = vcat(range(-U/2,U/2,length=N),range(0.1,1.0,length=N))
    fit = curve_fit(model_ED_real, ν_i, vcat(real(Δ_i),imag(Δ_i)), p0)
    ϵp = round.(fit.param[1:N], digits=15)
    Vp = round.(fit.param[(N+1):end], digits=15)
    check = abs(sum(Vp .^ 2) - 0.25)
    if check < best_check && minimum(abs.(Vp)) > Vp_min
        println("=> New best")
        global best_N = Nν_i
        global best_check = check
        global best_ϵp = ϵp
        global best_Vp = Vp
        global best_range = νrange 
        global best_fit = fit
        println("N = $Nν_i: sum(V^2_k) - 0.25 = $(round(check,digits=6))")
        println("ϵₗ = $(round.(best_ϵp, digits=4))")
        println("Vₗ = $(round.(best_Vp, digits=4))")

        println("    Anderson Parameter Checks: ")
        println("   ============================   ")
        println("1. min(|V_k|)       = $(round(minimum(abs.(best_Vp)), digits=4))")
        println("2. sum(V^2_k)       = $(round(sum(best_Vp .^ 2), digits=4))")
        println("3. min(|e_k|)       = $(round(minimum(abs.(best_ϵp)), digits=4))")
        println("4. min(|e_i - e_j|) = $(round(minimum(abs.(best_ϵp .- best_ϵp')  + Inf .* I ), digits=4))")
println("   ============================   ")
    end
end

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
