using thesis, PGFPlotsX
setup_pgfplotsx()


ts = LinRange(thesis.example_trange..., 21)
tts = LinRange(thesis.example_trange..., 201)

## Digital neuron
num_bits = 4
thresholds_dig = LinRange(thesis.example_yrange..., 2^num_bits+1)
levels_dig = 0.5*(thresholds_dig[1:end-1]+thresholds_dig[2:end])
ys_dig = thesis.example_input.(ts)
ys_bits_dig = thesis.discretize.(ys_dig, Ref(thresholds_dig), Ref(digits.(eachindex(thresholds_dig).-1, base=2, pad=num_bits)))
ys_num_dig = thesis.discretize.(ys_dig, Ref(thresholds_dig), Ref(levels_dig))
ys_trans_dig = thesis.discretize.(ys_dig, Ref(levels_dig), Ref(string.(0:2^num_bits,base=2^num_bits)))
y_rec_dig(t) = ys_num_dig[searchsortedlast(ts, t)]
τ = minimum(diff(ts))


β = 2.0
δ = 0.3
## LIF neuron
sol_lif,spikes_lif=lif_pair_get_solution(thesis.example_trange, thesis.example_input;α=β,β=β,δ=δ)
y_lif=thesis.example_input_integral.(spikes_lif)
y_rec_lif = t->sol_lif(t,idxs=2)
y_filtered = t->sol_lif(t,idxs=3)

δ = 0.15
## LNP neuron
sol_lnp,spikes_lnp=lnp_pair_get_solution(thesis.example_trange, thesis.example_input ;β=β,δ=δ)
y_lnp=thesis.example_input_integral.(spikes_lnp)
y_rec_lnp = t->sol_lnp(t,idxs=1)
y_filtered = t->sol_lnp(t,idxs=2)


## Build the monstrosity

@pgf gp = TikzPicture(
    GroupPlot(
    {
        group_style={group_size="3 by 4"},
        major_tick_style={draw="none"}, 
    },
    # input for digital
    {title="digital encoding", height="4cm",width="6cm", grid = "minor", ylabel=raw"input $s(t)$", xticklabel="\\empty", minor_xtick=ts, minor_ytick=levels_dig, ymin=thesis.example_yrange[1], ymax=thesis.example_yrange[2], xmin=thesis.example_trange[1], xmax=thesis.example_trange[2]},
    PlotInc(
        {no_marks, very_thick},
        Table(tts, thesis.example_input.(tts))
    ),
    (Plot(
        {no_markers, color="Set2-B", very_thick},
        Coordinates([(t,0),(t,thesis.example_input(t))])
        ) for (i,t) ∈ enumerate(ts))...,
    PlotInc(
        {only_marks, color="Set2-B"},
        Table(ts, ys_dig)
    ),
    PlotInc(
        {only_marks, mark="+", color="black"},
        Table(ts, ys_num_dig)
    ),
    [raw"\coordinate (bot11) at (axis cs:5,0);"],
    # input for lif
    {title="lIF encoding", height="4cm",width="6cm", grid = "minor", yminorgrids="false",xticklabel="\\empty", minor_xtick=spikes_lif, ymin=thesis.example_yrange[1], ymax=thesis.example_yrange[2], xmin=thesis.example_trange[1], xmax=thesis.example_trange[2]},
    PlotInc(
        {no_marks, very_thick, fill="Set2-B", fill_opacity=0.5},
        Table([thesis.example_trange[1];tts;thesis.example_trange[2]], [thesis.example_yrange[1];thesis.example_input.(tts);thesis.example_yrange[1]])
    ),
    # (Plot(
    #     {no_markers, color="black", very_thick},
    #     Coordinates([(t,0),(t,thesis.example_input(t))])
    #     ) for (i,t) ∈ enumerate(ts))...,
    [raw"\coordinate (bot12) at (axis cs:5,0);"],
    # input for lnp
    {title="LNP encoding", height="4cm",width="6cm", grid = "minor", yminorgrids="false",xticklabel="\\empty", minor_xtick=spikes_lnp, ymin=thesis.example_yrange[1], ymax=thesis.example_yrange[2], xmin=thesis.example_trange[1], xmax=thesis.example_trange[2]},
    PlotInc(
        {no_marks, very_thick, fill="Set2-B", fill_opacity=0.5},
        Table([thesis.example_trange[1];tts;thesis.example_trange[2]], [thesis.example_yrange[1];thesis.example_input.(tts);thesis.example_yrange[1]])
    ),
    # (Plot(
    #     {no_markers, color="black", very_thick},
    #     Coordinates([(t,0),(t,thesis.example_input(t))])
    #     ) for (i,t) ∈ enumerate(ts))...,
    [raw"\coordinate (bot13) at (axis cs:5,0);"],
    # mechanism for digital
    {height="4cm",width="6cm", xticklabel="\\empty", "group/empty plot"},
    # mechanism for lif
    {xtick_pos="bottom", yticklabel_pos="left", ytick_pos="right", tick_align="outside", color="black", minor_tick_length="2mm", height="4cm",width="6cm", grid = "minor", ylabel=raw"{$\int_{0}^t s(t)dt$}", xticklabel="\\empty", minor_xtick=spikes_lif, minor_ytick=y_lif, ymin=0, xmin=thesis.example_trange[1], xmax=thesis.example_trange[2]},
    PlotInc(
        {no_marks, very_thick, color="Set2-B"},
        Table(tts, thesis.example_input_integral.(tts))
    ),
    Plot(
        {only_marks, mark="+", color="black"},
        Table(spikes_lif, y_lif)
    ),
    [raw"\coordinate (top22) at (axis cs:5,\pgfkeysvalueof{/pgfplots/ymax});", raw"\coordinate (bot22) at (axis cs:5,0);"],
    # mechanism for lnp
    {xtick_pos="bottom", yticklabel_pos="left", ytick_pos="right", tick_align="outside", minor_tick_length="2mm", height="4cm",width="6cm", grid = "minor", xticklabel="\\empty", minor_xtick=spikes_lnp, minor_ytick=y_lnp, ymin=0, xmin=thesis.example_trange[1], xmax=thesis.example_trange[2]},
    PlotInc(
        {no_marks, very_thick, color="Set2-B"},
        Table(tts, thesis.example_input_integral.(tts))
    ),
    Plot(
        {only_marks, mark="+", color="black"},
        Table(spikes_lnp, y_lnp)
    ),
    [raw"\coordinate (top23) at (axis cs:5,\pgfkeysvalueof{/pgfplots/ymax});", raw"\coordinate (bot23) at (axis cs:5,0);"],
    # message for digital
    {height="2cm",width="6cm", grid = "minor", ylabel="message", xticklabel="\\empty", minor_xtick=ts, ytick="\\empty", ymin=0, ymax=1, xmin=thesis.example_trange[1], xmax=thesis.example_trange[2]},
    PlotInc(
        Table([],[])
    ),
    ["\\node[] at (axis cs:$t,0.5) {\\strut $y};" for (t,y) ∈ zip(0.5*(ts[1:end-1].+ts[2:end]), ys_trans_dig)],
    ["\\coordinate (label31) at (\$(axis cs:10,0)+(0,-4pt)\$);", raw"\coordinate (top31) at (axis cs:5,1);\coordinate (bot31) at (axis cs:5,0);"],
    # message for lif
    {height="2cm",width="6cm", grid = "minor", ylabel="spikes", xticklabel="\\empty", xtick="\\empty", ytick="\\empty", ymin=0, ymax=1, xmin=thesis.example_trange[1], xmax=thesis.example_trange[2]},
    (
        Plot(
            {no_markers, color="black", very_thick},
            Coordinates([(t,0),(t,1)])
        ) for (i,t) ∈ enumerate(spikes_lif)
    )...,
    ["\\coordinate (label32) at (\$(axis cs:10,0)+(0,-4pt)\$);", raw"\coordinate (top32) at (axis cs:5,\pgfkeysvalueof{/pgfplots/ymax});", raw"\coordinate (bot32) at (axis cs:5,0);"],
    # message for lnp
    {height="2cm",width="6cm", grid = "minor", xticklabel="\\empty", xtick="\\empty", ytick="\\empty", ymin=0, ymax=1, xmin=thesis.example_trange[1], xmax=thesis.example_trange[2]},
    (
        Plot(
            {no_markers, color="black", very_thick},
            Coordinates([(t,0),(t,1)])
        ) for (i,t) ∈ enumerate(spikes_lnp)
    )...,
    ["\\coordinate (label33) at (\$(axis cs:10,0)+(0,-4pt)\$);", raw"\coordinate (top33) at (axis cs:5,\pgfkeysvalueof{/pgfplots/ymax});", raw"\coordinate (bot33) at (axis cs:5,0);"],
    # reconstruction for digital
    {height="4cm",width="6cm", grid = "minor", xlabel=raw"time $t$", ylabel="reconstruction", minor_xtick=ts, minor_ytick=levels_dig, ymin=thesis.example_yrange[1], ymax=thesis.example_yrange[2], xmin=thesis.example_trange[1], xmax=thesis.example_trange[2]},
    PlotInc(
        {no_marks, very_thick},
        Table(tts, thesis.example_input.(tts))
    ),    
    Plot(
        {no_marks, very_thick, color="black"},
        Table(tts, y_rec_dig.(tts))
    ),
    PlotInc(
        {no_marks, color="Set2-B", very_thick},
        Table(tts, (thesis.example_input_integral.(tts).-thesis.example_input_integral.(tts.-τ))/τ)
    ),
    [raw"\coordinate (top41) at (axis cs:5,\pgfkeysvalueof{/pgfplots/ymax});", raw"\coordinate (bot42) at (axis cs:5,0);"],
    # reconstruction for lif
    {height="4cm",width="6cm", grid = "minor", yminorgrids="false",xlabel=raw"time $t$", minor_xtick=spikes_lif, ymin=thesis.example_yrange[1], ymax=thesis.example_yrange[2], xmin=thesis.example_trange[1], xmax=thesis.example_trange[2]},
    PlotInc(
        {no_marks, very_thick},
        Table(tts, thesis.example_input.(tts))
    ),  
    Plot(
        {no_marks, very_thick, color="black"},
        Table(tts, y_rec_lif.(tts))
    ),  
    PlotInc(
        {no_marks, color="Set2-B", very_thick},
        Table(tts, y_filtered.(tts))
    ),
    [raw"\coordinate (top42) at (axis cs:5,\pgfkeysvalueof{/pgfplots/ymax});", raw"\coordinate (bot42) at (axis cs:5,0);"],
    # reconstruction for lnp
    {height="4cm",width="6cm", grid = "minor", yminorgrids="false",xlabel=raw"time $t$", minor_xtick=spikes_lnp, ymin=thesis.example_yrange[1], ymax=thesis.example_yrange[2], xmin=thesis.example_trange[1], xmax=thesis.example_trange[2]},
    PlotInc(
        {no_marks, very_thick},
        Table(tts, thesis.example_input.(tts))
    ),    
    Plot(
        {no_marks, very_thick, color="black"},
        Table(tts, y_rec_lnp.(tts))
    ),
    PlotInc(
        {no_marks, color="Set2-B", very_thick},
        Table(tts, y_filtered.(tts))
    ),
    [raw"\coordinate (top43) at (axis cs:5,\pgfkeysvalueof{/pgfplots/ymax});", raw"\coordinate (bot43) at (axis cs:5,0);"],
    ), 
    raw"\draw[->, shorten >= 4pt, shorten <= 4pt] (bot11) -- node[draw, fill=white, rounded corners, align=center] {analog-to-digital\\conversion} (top31);",
    raw"\draw[->, shorten >= 4pt, shorten <= 4pt] (bot31) -- (top41);",
    raw"\draw[->, shorten >= 4pt, shorten <= 4pt] (bot12) -- (top22);",
    raw"\draw[->, shorten >= 4pt, shorten <= 4pt] (bot22) -- (top32);",
    raw"\draw[->, shorten >= 4pt, shorten <= 4pt] (bot32) -- (top42);",
    raw"\draw[->, shorten >= 4pt, shorten <= 4pt] (bot13) -- (top23);",
    raw"\draw[->, shorten >= 4pt, shorten <= 4pt] (bot23) -- (top33);",
    raw"\draw[->, shorten >= 4pt, shorten <= 4pt] (bot33) -- (top43);",
    "\\node[anchor=north east] at (label31) {($(sum(sum.(ys_bits_dig))) bits)};",
    "\\node[anchor=north east] at (label32) {($(length(spikes_lif)) spikes)};",
    "\\node[anchor=north east] at (label33) {($(length(spikes_lnp)) spikes)};",
)
#pgfsave("fig/encoding_schemes.pdf",gp)
