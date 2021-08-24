using thesis, PGFPlotsX, DifferentialEquations, DifferentialEquations.EnsembleAnalysis, Statistics, ProgressMeter
setup_pgfplotsx()

## Helper functions
"""the expected firing rate of an lIF neuron with threshold `δ` and time-constant `α` in response to a constant input `x`"""
r_est_lif(x,δ,α)::Float64 = x ≤ α*δ ? 0.0 : -α/log(1-δ*α/x)

"""the expected firing rate of an LNP neuron with gain `λ` and time-constant `α` in response to a constant input `x`"""
r_est_lnp(x,λ,α)::Float64 = x ≤ 0.0 ? 0.0 : λ*x

"""the expected membrane potential of the decoding lIF neuron for an encoding lIF neuron with threshold `δ` and time-constant `α` in response to a constant input `x`"""
z_est_lif(x,δ,α) = δ*r_est_lif(x,δ,α)

"""the expected threshold of an lIF neuron with time-constant `α` that leads to a firing rate `r` in response to a constant input `x`"""
δ_est(x,r, α) = x/α*(1-exp(-α/r))

"""the expected gain of an LNP neuron with time-constant `α` that leads to a firing rate `r` in response to a constant input `x`"""
λ_est(x,r, α) = r/x

"""the RMSE expected for a ternary value in a finite `range` """
RMSE_ternary(range) = √(1/12)*range/3

"""the RMSE expected for a binary-encoded value in a finite `range` with fixed number of `bits` """
RMSE_digital(bits, range) = √(1/12)*range*2.0^(-bits)

"""the RMSE expected for an lIF neuron with threshold δ and time-constant α """
RMSE_lif(δ, α) = √(1/12)*α*δ

"""the threshold δ expected for an lIF neuron with time-constant α required for a fixed `RMSE`"""
RMSE_lif_inverse(RMSE, α) = √(12)*RMSE/α

"""the RMSE expected for an LNP neuron with gain `λ` and time-constant α """
RMSE_lnp(λ, α) = √(α/(2*λ))

"""the gain `λ` expected for an LNP neuron with time-constant α required for a fixed `RMSE`"""
RMSE_lnp_inverse(RMSE, α) = α/(2*RMSE^2)

"""the estimated rate at which the expected RMSE will equal that of an `nbits` encoding """
r_nbit_lif(s₀, nbits, α) = r_est_lif(s₀, RMSE_lif_inverse(RMSE_digital(nbits, 1.0), α), α)

"""the estimated rate at which the expected RMSE will equal that of a ternary encoding """
r_ternary_lif(α) = r_est_lif(s₀, RMSE_lif_inverse(RMSE_ternary(1.0), α), α)

## Configuration types

abstract type SNNConfig end;

struct lIFConfig{A,D,II,W,IW} <: SNNConfig
    num_neurons::Int
    α::A
    δ::D
    input_neurons::II
    input_weights::IW
    weights::W
    can_spike::Vector{Bool}
end

struct LNPConfig{A,L,II,W,IW} <: SNNConfig
    num_neurons::Int
    α::A
    λ₀::L
    input_neurons::II
    input_weights::IW
    weights::W
    can_spike::Vector{Bool}
end

## ODE and problem definition for lIF and LNP neurons

"""dynamics of leaky neurons (lIF and LNP) with exponential filter"""
function f′(du, u, p, t)
    du .= -p.cfg.α*u

    # work-around for jump problems (broadcasting doesnt work)
    tmp = p.cfg.input_weights .* p.input(t)
    for idx ∈ p.cfg.input_neurons
        du[idx] += tmp[idx]
    end
end

function lIF_spike_condition(out, u, t, integrator)
    out .= u .- integrator.p.cfg.δ
    # ignore neurons that cant spike
    out[(~).(can_spike)] .= Inf
end

function lIF_spike!(integrator, idx)
    integrator.u .+= integrator.p.cfg.weights[:,idx]
    
    has_spiked = zeros(Bool, integrator.p.cfg.num_neurons)
    has_spiked[idx] = true
    push!(integrator.p.spike_train, (integrator.t, idx))
    while true
        new_spikes = findall(integrator.p.cfg.can_spike .* (~).(has_spiked) .* (integrator.u .>= integrator.p.cfg.δ))
        if isempty(new_spikes)
            break
        end
        for i ∈ new_spikes
            integrator.u .+= integrator.p.cfg.weights[:,i]
            has_spiked[i] = true
            push!(integrator.p.spike_train, (integrator.t, i))
        end
    end
    integrator.u[has_spiked] .= 0.0
end

function make_LNP_jump(i)
    function LNP_rate(u,p,t)
        return u[i]*p.cfg.λ₀[i]
    end

    function LNP_spike!(integrator)
        integrator.u.u .+= integrator.p.cfg.weights[:,i]
        push!(integrator.p.spike_train, (integrator.t, i))
    end

    VariableRateJump(LNP_rate, LNP_spike!)
end

"""run an ensemble of lIF encoder-decoder-pair simulations over timespan `tspan` for a given `input` (function of time) and an lIFConfig `p`
Returns a solution object that tracks, for each trajectory, the decoding neuron's membrane potential (first trace) and spike rate (second trace)
"""
function runEnsemble(tspan, input, p::lIFConfig; trajectories=100, kwargs...)
    p.weights[2,1]=p.δ[1] * p.α

    spike_trains = [Tuple{Float64,Int}[] for i ∈ 1:trajectories]

    prob_func(prob,i,repeat) = ODEProblem(f′, rand(p.num_neurons).*p.δ, tspan, (spike_train=spike_trains[i], cfg=p, input=input))
    prob = prob_func(nothing,1,0)
    cb = VectorContinuousCallback(lIF_spike_condition, lIF_spike!, nothing, p.num_neurons)
    ensemble_prob = EnsembleProblem(
        prob,
        prob_func=prob_func; 
        output_func=(sol,i)->(
            [
                sol(tspan[2])[2], 
                length(filter(((t,n),)->(t>(tspan[2]+tspan[1])/2), spike_trains[i]))
            ],
            false
        )
    )
    sol = solve(ensemble_prob, Tsit5(), scheduler(); callback=cb, trajectories=trajectories, kwargs...)

    return hcat(sol...)'
end

"""run an ensemble of LNP encoder-decoder-pair simulations over timespan `tspan` for a given `input` (function of time) and an LNPConfig `p`
Returns a solution object that tracks, for each trajectory, the decoding neuron's membrane potential (first trace) and spike rate (second trace)
"""
function runEnsemble(tspan, input, p::LNPConfig; trajectories=100, kwargs...)
    p.weights[2,1]=p.α / p.λ₀[1]

    spike_trains = [Tuple{Float64,Int}[] for i ∈ 1:trajectories]

    prob_func(prob,i,repeat) = JumpProblem(
        ODEProblem(f′, rand(p.num_neurons), tspan, (spike_train=spike_trains[i], cfg=p, input=input)), 
        Direct(),
        JumpSet((make_LNP_jump(j) for j ∈ 1:p.num_neurons if p.can_spike[j])...)
    )
    prob = prob_func(nothing,1,0)
    ensemble_prob = EnsembleProblem(
        prob,
        prob_func=prob_func; 
        output_func=(sol,i)->(
            [
                sol(tspan[2])[2], 
                length(filter(((t,n),)->(t>(tspan[2]+tspan[1])/2), spike_trains[i]))
            ],
            false
        )
    )
    sol = solve(ensemble_prob, Tsit5(), scheduler(); trajectories=trajectories, kwargs...)

    return hcat(sol...)'
end

"""simulates an ensemble of encoder-decoder pairs of neurons (lIF and LNP) for constant inputs with different `signal_levels`.
returns the `signal_levels`, means `μ` of lIF and LNP neurons and RMSE `err` of lIF and LNP neurons
"""
function calculate_transfer(tspan, signal_levels, lIF_cfg, LNP_cfg; kwargs...)
    μ = (lif=Float64[],lnp=Float64[])
    err = (lif=(neg=Float64[], pos=Float64[]),lnp=(neg=Float64[], pos=Float64[]))
    @showprogress 1.0 "Calculating transfer functions..." for level ∈ signal_levels
        for (i,cfg) ∈ enumerate((lIF_cfg, LNP_cfg))
            sol=runEnsemble(tspan, t->[level], cfg; kwargs...)
            
            quant = quantile(sol[:,1],[0.25,0.75])
            m = mean(sol[:,1])

            push!(err[i].neg,m-quant[1])
            push!(μ[i],m)
            push!(err[i].pos,quant[2]-m)
        end
    end
    return (signal_levels=signal_levels, μ=μ, err=err)
end


"""simulates an ensemble of encoder-decoder pairs of neurons (lIF and LNP) as well as a digital encoding for an `input` function.
returns digital encoding as well as decodings of lIF and LNP neurons 
"""
function compare_digital_lif_lnp(times, input, (bits, rate), lIF_cfg, LNP_cfg; kwargs...)
    tspan = (minimum(times),maximum(times))
    spike_trains = (digital=Tuple{Float64,Int}[], lif=Tuple{Float64,Int}[], lnp=Tuple{Float64,Int}[])
    sample_points = minimum(times):1.0/rate:maximum(times)
    samples = input.(sample_points)

    digital = map(t->round(samples[searchsortedlast(sample_points, t)]*2^bits)/2^bits, times)
    
    prob = ODEProblem(f′, rand(lIF_cfg.num_neurons).*lIF_cfg.δ, tspan, (spike_train=spike_trains.lif, cfg=lIF_cfg, input=input))
    cb = VectorContinuousCallback(lIF_spike_condition, lIF_spike!, nothing, lIF_cfg.num_neurons)
    sol_lIF = solve(prob, Tsit5(); callback=cb, kwargs...)
    
    prob = JumpProblem(
        ODEProblem(f′, rand(LNP_cfg.num_neurons), tspan, (spike_train=spike_trains.lnp, cfg=LNP_cfg, input=input)), 
        Direct(),
        JumpSet((make_LNP_jump(i) for i ∈ 1:LNP_cfg.num_neurons if LNP_cfg.can_spike[i])...)
    )
    sol_LNP = solve(prob, Tsit5(); kwargs...)

    return (solution=(digital=digital, lif=sol_lIF(times, idxs=2)[:].+0.5*lIF_cfg.δ[1]*lIF_cfg.α, lnp=sol_LNP(times, idxs=2)[:]), spikes=spike_trains)
end

"""simulates an ensemble of encoder-decoder pairs of neurons (lIF and LNP) for a fixed `input` function with different `scales` of firing rates (adjust δ or λ accordingly).
returns the `scales`, mean firing rates `rate` of lIF and LNP neurons and RMSE `rmse` of lIF and LNP neurons
"""
function calculate_error_vs_spikes!(tspan, scales, input, lIF_cfg, LNP_cfg; kwargs...)
    rate = (lif=Float64[],lnp=Float64[])
    rmse = (lif=Float64[],lnp=Float64[])
    @showprogress 1.0 "Calculating error/rate relation..." for s ∈ scales
        for (i,cfg) ∈ enumerate((lIF_cfg,LNP_cfg))
            if i==1
                cfg.δ[1] = δ_est(input(0)[1],s,cfg.α)
            else
                cfg.λ₀[1] = λ_est(input(0)[1],s,cfg.α)
            end
    
            sol=runEnsemble(tspan, input, cfg; kwargs...)
            
            push!(rate[i], mean(sol[:,2]./(tspan[2]-tspan[1])*2))
    
            v = sol[:,1]
            
            if i==1
                v.+= cfg.α.*cfg.δ[1]./2
            end
            push!(rmse[i], √(mean((v .- input(tspan[2])).^2)))
        end
    end

    return (scales=scales, rate=rate, rmse=rmse)
end

## Set up
use_multithreading = false                           # whether to use multithreading or not
num_neurons = 2                                     # we always simulate 2 neurons here
α=1000.0/25.0                                       # α = 40 corresponds to a time-scale of 25ms
δ=[δ_est(0.5,r_nbit_lif(0.5,4,α),α), 2.0]               # δ = 0.0015625 gives similar RMSE at a signal value 0.5 as a 4-bit encoding would
λ₀=1.0 ./ δ                                         # λ = 1/δ results in the same mean firing rate as for the lIF neuron
can_spike = [true,false]                            # only the encoding neuron can spike
input_neurons=[1]                                   # only the encoding neuron receives input
input_weights_lIF=[1.0]                             # the input-weights for the lIF neuron are 1
input_weights_LNP=input_weights_lIF.*α              # the input-weights for the LNP neuron are scaled up by α
weights=[0 0; 1 0].*δ'.*α                           # the weight from encoder to decode is scaled by αδ
trajectories = 50                                  # we simulate an ensemble of multiple trajectories to determine the distribution
signal_levels = LinRange(-0.2,1.0, 100)             # we evaluate the input-output function of the neurons for constant input signals with various levels
scales = 10.0 .^ LinRange(0, 4, 10)                 # we evaluate the RMSE for different parameterizations by sweeping over multiple scales (changes either δ or λ, respectively)
s₀ = 0.5                                            # we use a constant signal s(t)=0.5 to evaluate different parameterizations
equivalent_sampling_rate = α/2                      # a rectangular filter with equivalent delay as the exponential filter with time-constant α has width 2/α -> rate α/2
reference_bits = [2,4,6,8]                          # plot refrence lines for digital encodings using different numbers of bits
tspan = (0.0,10.0)                                  # we simulate for 10 seconds and take the last half of that for calculating averages
sine(t) = 0.4*sin(0.8*π*t)+0.5                      # we visualize the encoding for a sine-wave
plot_times = LinRange(-1,2.5,1000)                  # we plot the sine wave in the interval 0 to 10

inset_xmin,inset_xmax = 0.3,0.5
inset_ymin,inset_ymax = 0.7,1.0
plot_times_zoom = LinRange(inset_xmin,inset_xmax,200)

# the decoded signal will approximate the filtered version of the input signal 
#∫ᵗsin(ωτ)α exp(-α(t-τ))dτ 
#   = α exp(-αt)/(2i)∫ᵗ(exp((α+iω)τ)-exp((α-iω)τ))dτ 
#   = α exp(-αt)/(2i)[1/(α+iω)exp((α+iω)t)-1/(α-iω)exp((α-iω)t)]
#   = α 1/(2iα-2ω)exp(iωt)-1/(2iα+2ω)exp(-iωt)
#   = α (α sin(ωt) - ω cos(ωt))/(α² + ω²)
sine_filtered(t) = 0.4*(α*sin(0.8*π*t)-0.8*π*cos(0.8*π*t)) * α/(α^2+(0.8π)^2)+0.5              

##

scheduler = use_multithreading ? EnsembleThreads : EnsembleSerial
lIF_cfg = lIFConfig(num_neurons, α, copy(δ), input_neurons, copy(input_weights_lIF), copy(weights), can_spike)
LNP_cfg = LNPConfig(num_neurons, α, copy(λ₀), input_neurons, copy(input_weights_LNP), copy(weights), can_spike)

# Calculate the input-output transfer function
(s, μ, err) = calculate_transfer(tspan, signal_levels, lIF_cfg, LNP_cfg; trajectories=trajectories)

# Calculate input and decoding of a sine-wave
(decoded,spikes) = compare_digital_lif_lnp(plot_times, sine, (4, equivalent_sampling_rate), lIF_cfg, LNP_cfg; trajectories=trajectories)

# Calculate the RMSE/firing rate relationship as a function of varying δ or λ
(scales, rate, rmse) = calculate_error_vs_spikes!(tspan, scales, t->[s₀], lIF_cfg, LNP_cfg; trajectories=trajectories)

ss = 10 .^ LinRange(log10(minimum(scales)), log10(maximum(scales)), 200)

reference_rmse = [RMSE_digital(bits, 1.0) for bits in reference_bits]
reference_rate = [r_nbit_lif(s₀, bits,α) for bits in reference_bits]
reference_text = ["\\approx $(bits)-bit @ \$$(equivalent_sampling_rate)\$Hz " for bits in reference_bits]

t0=searchsortedfirst(plot_times,0)

## Plot the monstrosity
@pgf gp = TikzPicture(
    {remember_picture},  
    GroupPlot(
    {
        group_style={
            group_size="2 by 3",
            horizontal_sep="2cm",
            vertical_sep="2cm",
        },
        height="6cm",
        width="8cm",
    },
    {xlabel=raw"input value $c$", ylabel=raw"rate$(c)$", title=raw"decoded rate", legend_pos="north west", xmin=signal_levels[1], xmax=signal_levels[end]},
    Plot(
        {"name path=lnplow", draw="none", forget_plot},
        Table(signal_levels, μ.lnp .- err.lnp.neg)
    ),
    Plot(
        {"name path=lnphigh", draw="none", forget_plot},
        Table(signal_levels, μ.lnp .+ err.lnp.pos)
    ),
    PlotInc(
        {color="Set1-B", opacity=0.5, forget_plot},
        raw"fill between [of=lnplow and lnphigh]"
    ),
    PlotInc(
        {color="Set1-B", very_thick},
        Table(signal_levels, μ.lnp)
    ),
    LegendEntry("LNP response"),
    Plot(
        {"name path=liflow", draw="none", forget_plot},
        Table(signal_levels, μ.lif .- err.lif.neg)
    ),
    Plot(
        {"name path=lifhigh", draw="none", forget_plot},
        Table(signal_levels, μ.lif .+ err.lif.pos)
    ),
    PlotInc(
        {color="Set1-A", forget_plot},
        raw"fill between [of=liflow and lifhigh]"
    ),
    PlotInc(
        {color="Set1-A", very_thick},
        Table(signal_levels, μ.lif)
    ),
    LegendEntry("lIF response"),
    PlotInc(
        {dashed, very_thick},
        Table(signal_levels, max.(0,signal_levels.- 0.5*α*δ[1]))
    ),
    LegendEntry(raw"shifted ReLU"),
    PlotInc(
        {dotted, very_thick},
        Table(signal_levels, z_est_lif.(signal_levels, δ[1], α))
    ),
    LegendEntry(raw"expected lIF response"),
    {xlabel=raw"firing rate $r$ [Hz]", ylabel=raw"RMSE$(r)$", title="RMSE for constant input \$x(t) = $(s₀)\$", ymode="log", xmode="log",xmin=ss[1],xmax=ss[end]},
    (Plot(
        {forget_plot, thick, color= i==2 ? "Set1-C" : "black!30"},
        Coordinates([(ss[1],e),(ss[end],e)])
    ) for (i,e) ∈ enumerate(reference_rmse))...,
    PlotInc(
        {forget_plot, very_thick, color="Set1-B"},
        Table(ss, RMSE_lnp.(λ_est.(s₀,ss,α),α))
    ),
    PlotInc(
        {only_marks, very_thick, color="Set1-B"},
        Table(rate.lnp, rmse.lnp)
    ),
    LegendEntry(raw"LNP neuron"),
    PlotInc(
        {forget_plot, very_thick, color="Set1-A"},
        Table(ss, RMSE_lif.(δ_est.(s₀,ss,α),α))
    ),
    PlotInc(
        {only_marks, very_thick, color="Set1-A"},
        Table(rate.lif, rmse.lif)
    ),
    LegendEntry(raw"lIF neuron"),
    Plot(
        {only_marks, mark="+", color="black"},
        Table(reference_rate, reference_rmse)
    ),
    ["\\node[anchor=north east, color=$(i==2 ? "Set1-C" : "black!50")] at (axis cs:$(r),$(e)) {$(t)};" for (i,(r,e,t)) ∈ enumerate(zip(reference_rate, reference_rmse, reference_text))],
    {width="16.5cm", xshift="4.25cm", title="digital code vs. lIF neuron vs. LNP neuron", ylabel="decoded signal", xticklabels="\\empty", xmin=0,xmax=2.5,ymin=0, ymax=2.0},
    "\\coordinate (boxsw) at (axis cs:$(inset_xmin),$(inset_ymin));",
    "\\coordinate (boxse) at (axis cs:$(inset_xmax),$(inset_ymin));",
    "\\coordinate (boxne) at (axis cs:$(inset_xmax),$(inset_ymax));",
    "\\coordinate (boxnw) at (axis cs:$(inset_xmin),$(inset_ymax));",
    "\\filldraw[color=black!90, very thick, fill=black!10] (boxsw) rectangle (boxne);",
    PlotInc(
        {very_thick, color="Set1-B"},
        Table(plot_times,decoded.lnp)
    ),
    LegendEntry("LNP decoding ($(length(filter(>(0), first.(spikes.lnp)))) spikes)"),
    PlotInc(
        {very_thick, color="Set1-A"},
        Table(plot_times,decoded.lif)
    ),
    LegendEntry("lIF decoding ($(length(filter(>(0), first.(spikes.lif)))) spikes)"),
    PlotInc(
        {very_thick},
        Table(plot_times,decoded.digital)
    ),
    LegendEntry("digital decoding ($(sum(sum.(thesis.discretize.(decoded.digital, Ref(LinRange(0,1,2^4+1)), Ref(digits.(0:2^4,base=2,pad=4)))[t0:end])))) bits)"),
    Plot(
        {ultra_thick, dashed},
        Table(plot_times,sine_filtered.(plot_times))
    ),
    LegendEntry(raw"target signal $(s*\kappa_\alpha)(t)$"),
    raw"\coordinate (insetpos) at (axis cs:0.1,1.1);",
    {"group/empty plot", "axis background/.style={fill=none}", xmin=0,xmax=1,ymin=0,ymax=1},
    {yshift="1.75cm", width="16.5cm", height="3cm", xlabel="time [s]", xmin=0,xmax=2.5,ymin=0, ymax=2.0, ylabel="spikes", ytick="\\empty"},
    (Plot(
        {color="Set1-A", forget_plot, opacity=0.5},
        Coordinates([(t,0.05),(t,0.95)])
    ) for t ∈ first.(spikes.lif))...,
    (Plot(
        {color="Set1-B", forget_plot, opacity=0.5},
        Coordinates([(t,1.05),(t,1.95)])
    ) for t ∈ first.(spikes.lnp))...,
    ),
    Axis(
        {very_thick, height="1.75cm", width="5cm", scale_only_axis, "axis background/.style"={fill="black!10", opacity=0.9},  xmin=inset_xmin,xmax=inset_xmax, ymin=inset_ymin, ymax=inset_ymax, at={"(insetpos)"}, anchor={"south west"}, ticks="none"},
        PlotInc(
            {ultra_thick, color="Set1-B"},
            Table(plot_times,decoded.lnp)
        ),
        PlotInc(
            {ultra_thick, color="Set1-A"},
            Table(plot_times,decoded.lif)
        ),
        PlotInc(
            {ultra_thick},
            Table(plot_times,decoded.digital)
        ),
        Plot(
            {ultra_thick, dashed},
            Table(plot_times,sine_filtered.(plot_times))
        ),
        "\\coordinate (insw) at (axis cs:$(inset_xmin),$(inset_ymin));",
        "\\coordinate (inse) at (axis cs:$(inset_xmax),$(inset_ymin));",
        "\\coordinate (inne) at (axis cs:$(inset_xmax),$(inset_ymax));",
        "\\coordinate (innw) at (axis cs:$(inset_xmin),$(inset_ymax));",
    ),
    raw"\draw[color=black!75, very thick] (boxsw) -- (insw);",
    raw"\draw[color=black!75, very thick] (boxse) -- (inse);",
)
#pgfsave("fig/encoding_comparison.pdf",gp)
