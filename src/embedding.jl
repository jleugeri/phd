using DifferentialEquations, thesis

function f′(du, u, p, t)
    du .= -p.α .* u
    du[1] += p.input(t)
    du[2:end] .+= p.α*u[1:end-1]
end

p = (
    input=t-> 0.5*(sin(t)) + 0.3*(cos(0.7+1.5*t)),
    α = 2.0,
    δ = zeros(5),
    spiketimes = [Float64[] for i ∈ 1:5]
)
tspan = (0.0,10.0)
u₀ = copy(p.δ)
cb = CallbackSet((ContinuousCallback((u,t,integrator) -> u[i]-integrator.p.δ[i], integrator->push!(integrator.p.spiketimes[i], integrator.t)) for i ∈ eachindex(p.δ))...)

sol = solve(ODEProblem(f′, u₀, tspan.-(10.0,0), p); callback=cb)

p_in = plot(sol, tspan=tspan, xlims=tspan.+(0,0.2), xlabel=false, ylabel="signal embedding")
plot!(p.input, tspan..., color=:black, label="s(t)")
hline!([p.δ], color=thesis_colors[end], linestyle=:dash, label="threshold")

p_spike = plot(xlabel="time", ylabel="population code")
for (i,spiketrain) ∈ enumerate(p.spiketimes)
    spiketrain = filter(s->tspan[1]<=s<=tspan[2], spiketrain)
    scatter!(p_in, spiketrain, fill(p.δ[i],length(spiketrain)), markersize=3, markercolor=:black, markerstrokecolor=i, markerstrokewidth=2, label="")
    plot!(p_spike, [spiketrain'; spiketrain'], [fill(i,length(spiketrain))';fill(i+0.75,length(spiketrain))'], linewidth=2, color=i, label="")
    scatter!(p_spike, spiketrain, fill(i+0.75,length(spiketrain)), markersize=3,  markercolor=:black, markerstrokecolor=i, markerstrokewidth=2, label="")
end
plot!(p_spike, xlims=tspan.+(0,0.2), ylims=[1.0,6.25], yticks=1:5, label=false)

plot(p_in, p_spike, layout=grid(2,1, heights=[0.7,0.3]))

savefig("figures/embedding.pdf")
savefig("figures/embedding.svg")
