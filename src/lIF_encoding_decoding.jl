using Plots, DifferentialEquations
gr()

function f′(du, u, p, t)
    du[1] = p.input(t) - p.α*u[1]
    du[2] = -p.β * u[2]
    du[3] = p.β*(p.input(t) - u[3])
end

spike_condition(u,t,integrator) = u[1]-integrator.p.δ

function spike!(integrator)
    integrator.u[1] = 0.0
    integrator.u[2] += 1.0
    push!(integrator.p.spiketimes, integrator.t)
end

u₀ = [0.0,1.0,1.0]
tspan = (0.0,10.0)
p = (
    input=t-> 0.5*(sin(t)+1) + 0.3*(cos(0.7+1.5*t)+1),
    α = 0.5,
    β = 0.5,
    δ = 0.2,
    spiketimes = Float64[0.0]
)

cb = ContinuousCallback(spike_condition, spike!, nothing)
sol = solve(ODEProblem(f′, u₀, tspan, p); callback=cb)

function make_analytical_solution(spike_times; α=1.0, δ=0.2, basis=:sample_and_hold)
    return function analytical_reconstruction(t)
        i = searchsortedfirst(spike_times,t)
        return if i > length(spike_times)
            0.0
        else
            t₁ = spike_times[i]
            t₀ = i<=1 ? -Inf : spike_times[i-1]
            dt = t₁-t₀
            norm = √((1.0-exp(-2.0α*dt))/(2.0α))
            f = δ/norm # fourier coefficient

            if basis == :sample_and_hold
                1.0/√(dt) * f                
            else#if basis == :truncated_exponential
                exp(-α*(t₁-t))/norm * f
            end
        end
    end
end

ŝ = make_analytical_solution(p.spiketimes; α=p.α, δ=p.δ, basis=:sample_and_hold)
mean_dt = (tspan[2]-tspan[1])/length(p.spiketimes)
mean_scale = 2*p.δ*p.β/(1+exp(-p.α*mean_dt))

begin
p1=vline(p.spiketimes, linewidth=2, color="#babdb6ff", label="spikes ", xlims=(0,10.5), ylims=(0,1.7))
plot!(p.input, tspan..., linewidth=2, color="#729fcfff", label="input ")
plot!(t->sol(t,idxs=3), tspan..., linewidth=2, color="#fcaf3eff", linestyle=:dash, label="filtered input  ")
plot!(ŝ, tspan..., linewidth=2, color="#729fcfff", linestyle=:dash, label="acausal rec.  ")
plot!(t->sol(t,idxs=2)*mean_scale+exp(-p.β*t)*(1-mean_scale), tspan..., linewidth=2, color="#fcaf3eff", label="causal LTI rec.  ")
plot!(xlabel="time", legend=false)

idx = searchsortedlast(p.spiketimes, 5.0)
(t₀,t₁,t₂) = p.spiketimes[idx:idx+2]
xlims = (t₀-0.1, t₂+0.1)

p2 = vline(p.spiketimes, linewidth=2, color="#babdb6ff", xlims=xlims, ylims=(0,1.0), label="spikes ", legend=false)
plot!(p.input, xlims..., linewidth=2, color="#729fcfff", label="input ")
plot!(t->sol(t,idxs=3), tspan..., linewidth=2, color="#fcaf3eff", linestyle=:dash, label="filtered input ")
plot!(t->sol(t,idxs=1), xlims..., linewidth=2, color="#ad7fa8ff", label="state ")
hline!([p.δ], color="#ef2929ff", linewidth=2, linestyle=:dash, label="threshold ")
plot!(t->sol(t,idxs=2)*mean_scale+exp(-p.β*t)*(1-mean_scale), xlims..., linewidth=2, color="#fcaf3eff", label="causal LTI rec.")
plot!(t->sol(t,idxs=2)*mean_scale+exp(-p.β*t)*(1-mean_scale)-exp(-p.β*(t-t₁))*(t>t₁)*mean_scale, t₁, xlims[end], linestyle=:dot, linewidth=2, color="#fcaf3eff", label="")

p3 = plot(showaxis=false, framestyle=:none, showgrid=false)
plot!([],[], linewidth=2, color="#babdb6ff", label="spikes ")
plot!([],[], linewidth=2, color="#729fcfff", label="input ")
plot!([],[], linewidth=2, color="#fcaf3eff", linestyle=:dash, label="filtered input ")
plot!([],[], linewidth=2, color="#ad7fa8ff", label="state ")
plot!([],[], linewidth=2, color="#ef2929ff", linestyle=:dash, label="threshold ")
plot!([],[], linewidth=2, color="#729fcfff", linestyle=:dash, label="acausal rec.  ")
plot!([],[], linewidth=2, color="#fcaf3eff", label="causal LTI rec. ")
plot!(legend=(0.0, 0.9))

l = @layout [   a{0.7w} b
c{0.5h}     ]

plot(p2,p3,p1, layout=l)
end

savefig("figures/reconstruction_lif.pdf")
savefig("figures/reconstruction_lif.svg")