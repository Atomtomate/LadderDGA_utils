using LsqFit
using LinearAlgebra
using HDF5
using NLsolve
using FiniteDiff
using NLSolvers
using Dispersions

#f = h5open(ARGS[1])
N = 4#tryparse(Int,ARGS[2])
Nν_max = 1#tryparse(Int,ARGS[3])

β = read_attribute(f[".config"], "general.beta")
U = read_attribute(f[".config"], "atoms.1.udd")
μ = read(f["dmft-last/mu/value"])
Δ_in_tmp = read(f["dmft-last/ineq-001/fiw/value"])
G0_in_tmp = read(f["dmft-last/ineq-001/g0iw/value"]) 

Δ_in = (Δ_in_tmp[:,1,1] .+ Δ_in_tmp[:,2,1]) ./ 2
G0_in = (G0_in_tmp[:,1,1] .+ G0_in_tmp[:,2,1]) ./ 2

iν_in = 1im .* read(f[".axes/iw"])

nh = floor(Int, length(iν_in)/2 + 1)

kG = gen_kGrid("2Dsc-0.25-0.05-0.025", 100)

function model_ED_Δ(iν, p)
    Δ_fit = zeros(ComplexF64,length(iν))
    for (i,νn) in enumerate(iν)
        Δ_fit[i] = sum((p[(N+1):end] .^ 2) ./ (νn .- p[1:N]))
    end
    return conj.(Δ_fit)
end

function model_ED_G0(iν, p)
    Δ_fit = zeros(ComplexF64,length(iν))
    for (i,νn) in enumerate(iν)
        Δ_fit[i] = sum((p[(N+1):end] .^ 2) ./ (νn .- p[1:N]))
    end
    return [kintegrate(kG, 1 ./ (iν[i] .+ 0 .+ dispersion(kG) .- conj.(Δ_fit)[i])) for i in 1:length(iν)]
end


function model_ED_Δ_real(iν, p)
    Δ_fit = model_ED_Δ(iν,p)
    return vcat(real(Δ_fit),imag(Δ_fit))
end

function model_ED_G0_real(iν, p)
    G0_fit = model_ED_G0(iν,p)
    return vcat(real(G0_fit),imag(G0_fit))
end



best_N = 0
best_check = Inf
best_fit = nothing
best_ϵp = nothing
best_Vp = nothing
best_range = nothing
Vp_min = 0.010

println("Trying to determine best fit for NBath = $N, U = $U, β=$β, min(abs(Vp)) = $Vp_min in range of fermionic frequencies $(10:(10+Nν_max))")

for Nν_i in 20:(20+Nν_max)
    # Variable definitions
    νrange = (nh-Nν_i):(nh+Nν_i-1)
    Δ_i  = Δ_in[νrange]
    G0_i = G0_in[νrange]
    ν_i  = iν_in[νrange]

    # Objective function definitions
    tmp = similar(iν_in) 
    function objective(p::Vector; expN::Int = 1)
        GWeiss!(tmp, iν_in, μ, p[1:N], p[N+1:end])
        sum(abs.(tmp .- G0_i))/(abs.(ν_i) .^ expN)
    end
    function grad(∇f, x) 
        ∇f = FiniteDiff.finite_difference_gradient(objective,x)
        return ∇f
    end

    function hess(∇²f, x) 
        ∇²f = FiniteDiff.finite_difference_hessian(objective,x)
        return ∇²f
    end
    objective_grad(∇f,x) = objective(x), grad(∇f, x)
    function objective_grad_hess(∇f,∇²f,x) 
        f, ∇f = objective_grad(∇f,x) 
        ∇²f   = hess(∇²f, x)
        return f,∇f,∇²f
    end 
    
    # G0 

    # Least Squares G0
    p0_G0 = vcat(range(-U/2,U/2,length=N),range(0.1,1.0,length=N))
    fit_G0_lsq = curve_fit(model_ED_G0_real, ν_i, vcat(real(G0_i),imag(G0_i)), p0_G0)
    ϵp = round.(fit_G0_lsq.param[1:N], digits=15)
    Vp = round.(fit_G0_lsq.param[(N+1):end], digits=15)
    check = abs(sum(Vp .^ 2) - 0.25)
        println("=> CHECKS (G0 Least Squares)")
        println("N = $Nν_i: sum(V^2_k) - 0.25 = $(round(check,digits=6))")
        println("ϵₗ = $(round.(ϵp, digits=4))")
        println("Vₗ = $(round.(Vp, digits=4))")

        println("    Anderson Parameter Checks: ")
        println("   ============================   ")
        println("1. min(|V_k|)       = $(round(minimum(abs.(Vp)), digits=4))")
        println("2. sum(V^2_k)       = $(round(sum(Vp .^ 2), digits=4))")
        println("3. min(|e_k|)       = $(round(minimum(abs.(ϵp)), digits=4))")
        println("4. min(|e_i - e_j|) = $(round(minimum(abs.(ϵp .- ϵp')  + Inf .* I ), digits=4))")
println("   ============================   ")
#
    # Least Squares Δ
    p0_Δ = vcat(range(-U/2,U/2,length=N),range(0.1,1.0,length=N))
    fit_Δ_lsq = curve_fit(model_ED_Δ_real, ν_i, vcat(real(Δ_i),imag(Δ_i)), p0_Δ)
    ϵp = round.(fit_Δ_lsq.param[1:N], digits=15)
    Vp = round.(fit_Δ_lsq.param[(N+1):end], digits=15)
    check = abs(sum(Vp .^ 2) - 0.25)
        println("=> CHECKS (Δ Least Squares)")
        println("N = $Nν_i: sum(V^2_k) - 0.25 = $(round(check,digits=6))")
        println("ϵₗ = $(round.(ϵp, digits=4))")
        println("Vₗ = $(round.(Vp, digits=4))")

        println("    Anderson Parameter Checks: ")
        println("   ============================   ")
        println("1. min(|V_k|)       = $(round(minimum(abs.(Vp)), digits=4))")
        println("2. sum(V^2_k)       = $(round(sum(Vp .^ 2), digits=4))")
        println("3. min(|e_k|)       = $(round(minimum(abs.(ϵp)), digits=4))")
        println("4. min(|e_i - e_j|) = $(round(minimum(abs.(ϵp .- ϵp')  + Inf .* I ), digits=4))")
println("   ============================   ")
end

