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
    input=t-> 0.5*(sin(t)+1) + 0.3*(cos(0.7+1.5*t)+1)-0.3,
    α = 0.5,
    β = 0.5,
    δ = 0.2,
    spiketimes = Float64[0.0]
)

cb = ContinuousCallback(spike_condition, spike!, nothing)
sol = solve(ODEProblem(f′, u₀, tspan, p); callback=cb)


#ŝ = make_analytical_solution(p.spiketimes; α=p.α, δ=p.δ, basis=:sample_and_hold)
mean_dt = (tspan[2]-tspan[1])/length(p.spiketimes)
mean_scale = 2*p.δ*p.β/(1+exp(-p.α*mean_dt))

sample_times = LinRange(tspan..., length(p.spiketimes))
sample_val = p.input.(sample_times)
digitized_signal(t) = get(sample_val,searchsortedlast(sample_times,t),0.0)
q = (
    input=digitized_signal,
    α = 0.5,
    β = 0.5,
    δ = 0.2,
    spiketimes = Float64[0.0]
)
digital_sol = solve(ODEProblem(f′, u₀, tspan, q))

let last_spike=tspan[1], next_spike=tspan[1], spike_idx=0
    ts = LinRange(tspan...,300)[2:end]
    anim = @animate for (i,t) in enumerate(ts)
        if t >= next_spike
            last_spike = next_spike
            if spike_idx < length(p.spiketimes)
                spike_idx += 1
                next_spike = p.spiketimes[spike_idx]
            else
                next_spike = Inf
            end
        end

        sample_idx=searchsortedlast(sample_times,t)

        p_in_analog=plot(p.input, tspan[1],t, xlims=tspan, ylims=(-0.2,1.3), title="input", ylabel="analog")
        p_state_analog=plot(p.input, tspan[1],t, xlims=tspan, ylims=(-0.2,1.3), title="encoding")
        p_trans_analog=plot(p.input, tspan[1],t, xlims=tspan, ylims=(-0.2,1.3), title="transmission")
        p_out_analog=plot(t->sol(t,idxs=3), tspan[1],t, xlims=tspan, ylims=(-0.2,1.3), title="decoding")

        p_in_digital=plot(p.input, tspan[1],t, xlims=tspan, ylims=(-0.2,1.3), ylabel="digital")
        plot!(sample_times[1:sample_idx], p.input.(sample_times[1:sample_idx]), line=:stem, xlims=tspan, ylims=(-0.5,1.6))
        p_state_digital=plot(sample_times[1:sample_idx], p.input.(sample_times[1:sample_idx]), seriestype=:steppre, xlims=tspan, ylims=(-0.2,1.3))
        p_trans_digital=plot(sample_times[1:sample_idx], p.input.(sample_times[1:sample_idx]), marker=:o, line=:stem, xlims=tspan, ylims=(-0.2,1.3))
        t_samp = [sample_times[1:sample_idx]; t]
        val_samp= [sample_val[1:sample_idx];sample_val[sample_idx]]
        p_out_digital=plot(t->digital_sol(t,idxs=3), tspan[1], t, xlims=tspan, ylims=(-0.2,1.3))

        p_in_spike=plot(p.input, tspan[1],t, xlims=tspan, ylims=(-0.2,1.3), fill=(0,:gray), ylabel="spiking")
        plot!(p.input, last_spike,t, xlims=tspan, ylims=(-0.2,1.3), fill=(0,"#fcaf3eff"))
        p_state_spike=plot(t->sol(t,idxs=1), tspan[1],t, xlims=tspan, ylims=(-0.05  ,1.2*p.δ), color="#fcaf3eff")
        hline!([p.δ], color=:red, linestyle=:dash)
        p_trans_spike=vline(p.spiketimes[1:spike_idx], xlims=tspan)
        p_out_spike = plot(t->sol(t,idxs=2)*mean_scale+exp(-p.β*t)*(1-mean_scale), tspan[1],t, xlims=tspan, ylims=(-0.2,1.3))


        plot(p_in_analog, p_state_analog, p_trans_analog, p_out_analog,
             p_in_digital, p_state_digital, p_trans_digital, p_out_digital,
             p_in_spike, p_state_spike, p_trans_spike, p_out_spike, layout=grid(3,4), legend=false, grid=false, size=(1600,800), linewidth=2, fontsize=20)
    end
    gif(anim, "figures/lif_animation.gif", fps = 30)
end
