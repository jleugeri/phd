using thesis, Distributed, Plots, DifferentialEquations
thesis.init()

function f′(du, u, p, t)
    du[1] = -p.β * u[1]
end

function g′(du,u,p,t)
    du[1] = p.β*(p.input(t) - u[1])
    du[2] = 2p.β*(p.input(t)*p.β/(2p.λ₀) - u[2])
end

function spike!(integrator)
    integrator.u[1] += integrator.p.β/integrator.p.λ₀
#    push!(integrator.p.spiketimes, integrator.t)
end

function randomize_initial(prob, i, repeat)
    prob.prob.u0 .+= clamp.(randn(1).*sqrt(prob.prob.p.input(0.0) * prob.prob.p.β / (2*prob.prob.p.λ₀)), 0, Inf)
    prob
end

const tspan = (0.0,10.0)
const inp = t-> 0.5*(sin(t)+1) + 0.3*(cos(0.7+1.5*t)+1)
const λ₀s = [2.0, 10.0, 50.0]


function mkplot(i,λ₀)
    p = (
        input=inp,
        β = 0.5,
        spiketimes = Float64[0.0],
        λ₀= λ₀
    )

    u₀_mc = [p.input(0)]
    u₀_th=[p.input(0), p.input(0.0) * p.β / (2*p.λ₀)]


    prob_mc = MonteCarloProblem(JumpProblem(ODEProblem(f′, u₀_mc, tspan, p), Direct(), ConstantRateJump((u,p,t)->p.λ₀*p.input(t), spike!)), prob_func=randomize_initial)
    prob_th = ODEProblem(g′, u₀_th, tspan, p)

    xs = LinRange(tspan...,201)

    sol = solve(prob_mc, Tsit5(), saveat=xs)
    sol_theoretical = solve(prob_th, Tsit5(), saveat=xs)

    m,v=timeseries_point_meanvar(sol, xs)
    m̂,v̂=sol_theoretical(xs, idxs=1), sol_theoretical(xs, idxs=2)

    pp=plot(xs, m[1,:], ribbon=sqrt.(v[1,:]), color=:gray, title="λ=$(λ₀)")
    plot!(sol, idxs=1:10, linewidth=2, alpha=0.5, color=:gray)
    plot!(p.input, tspan..., color=1)
    plot!(xs, m̂ .+ sqrt.(v̂).*[1 -1], color=4, style=:dash, linewidth=2)
    plot!(m̂, color=4, linewidth=2, xlims=(0,10.5), ylims=(0,1.7), xlabel= i==3 ? "time" : "")
    
    if i!=3
        plot!(xticks=([],[]))
    end

    return pp
end

plts = Plots.Plot[]
for (i,λ₀) ∈ enumerate(λ₀s)
    push!(plts,mkplot(i,λ₀))
end

plot(plts..., layout=grid(3,1), size=(600,800))

savefig("fig/reconstruction_lnp.svg")
savefig("fig/reconstruction_lnp.pdf")