using DifferentialEquations, PGFPlotsX, LaTeXStrings, Colors
export lif_pair_get_solution, lnp_pair_get_solution, setup_pgfplotsx

setup_pgfplotsx() = begin
push!(PGFPlotsX.CUSTOM_PREAMBLE, 
raw"""
\usepackage{amsmath}
\usepgfplotslibrary{colorbrewer}
\usepgfplotslibrary{fillbetween}
\pgfplotsset{
    cycle list/Set1,
}
"""
)
end

setup_pgfplotsx()

example_trange=(0,10.0)
example_yrange=(0,2.0)
example_input=t-> 0.5*(sin(t)+1) + 0.3*(cos(0.7+1.5*t)+1)
example_input_integral=t-> -0.5*cos(t) + 0.3*sin(0.7+1.5*t)/1.5 + 0.8*t

discretize(x, levels, values) = values[searchsortedlast(levels, x)]

function lif_pair_f′(du, u, p, t)
    du[1] = p.α*p.input(t) - p.α*u[1]
    du[2] = -p.β * u[2]
    du[3] = p.β*(p.input(t) - u[3])
end

lif_pair_spike_condition(u,t,integrator) = u[1]-integrator.p.δ

function lif_pair_spike!(integrator)
    integrator.u[1] = 0.0
    integrator.u[2] += integrator.p.δ
    push!(integrator.p.spiketimes, integrator.t)
end

function lif_pair_get_solution(tspan, input;  α = 0.5, β = 0.5, δ = 0.2, u₀ = [0.0,1.0,1.0])
    spiketimes=Float64[]
    p = (;input,α,β,δ,spiketimes)
    cb = ContinuousCallback(lif_pair_spike_condition, lif_pair_spike!, nothing)
    sol = solve(ODEProblem(lif_pair_f′, u₀, tspan, p); callback=cb)
    return sol, spiketimes
end

function lnp_pair_f′(du, u, p, t)
    du[1] = -p.β * u[1]
    du[2] = p.β*(p.input(t) - u[2])
end

function lnp_pair_spike!(integrator)
    integrator.u[1] += integrator.p.δ
    push!(integrator.p.spiketimes, integrator.t)
end

function lnp_pair_get_solution(tspan, input;  β = 0.5, δ=0.2, u₀ = [1.0,1.0], kwargs...)
    spiketimes = Float64[]
    p = (;input,β, δ, spiketimes)
    rate_wrapped(u,p,t) = input(t)*β/δ
    prob = ODEProblem(lnp_pair_f′, u₀, tspan, p)
    jump = VariableRateJump(rate_wrapped, lnp_pair_spike!)
    jump_prob=JumpProblem(prob, Direct(),jump)
    sol = solve(jump_prob,Tsit5(); kwargs...)

    return sol, spiketimes
end
