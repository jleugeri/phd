using PGFPlotsX, Distributions, Cubature, LaTeXStrings, Interpolations

""" p(y|x) = N(f(x),σ)"""
p_y_x(y,x; σ=σ, slope=slope, context...) = pdf(Normal(f(x; slope=slope),σ),y)

""" p(x) = N(0,θ)"""
p_x(x; θ=θ, context...) = pdf(Normal(0,θ),x)

""" p(x,y) = p(y|x)p(x)"""
p_xy(x,y; context...) = p_y_x(y,x; context...)*p_x(x; context...)

""" p(y) = ∫ p(x,y) dx"""
p_y(y; θ=θ, context...) = hquadrature(x->p_xy(x,y; θ=θ, context...), -10*θ, 10*θ; reltol=ε)[1]

""" p(x|y) = p(x,y)/p(y)"""
p_x_y(x,y; context...) = p_xy(x,y; context...)/p_y(y; context...)

""" p(ŷ|y) = N(ŷ,σ)p(ŷ)/p(y)"""
p_ŷ_y(ŷ, y; context...) = p_x_y(f⁻(ŷ; context...),y; context...)*abs(f⁻′(ŷ; context...))

# How much does y reveal about x?
# Q: How much better than the prior do we know where x is after observing y?
# A: KLD from posterior to prior
# = ∫ p(y) ∫ p(x|y)log(p(x|y)/p(x)) dx dy
# = ∫∫ p(x,y) log(p(x|y)/p(x)) dxdy
# = ∫∫ p(x,y) log(p(x,y)/(p(x)p(y))) dxdy
# = MI(X,Y)

function observed_information(x,y; context...)
    P_x = p_x(x; context...)
    P_y = p_y(y; context...)
    P_y_x = p_y_x(y,x; context...)
    P_xy = P_y_x*P_x
    return P_xy*(log2(P_y_x)-log2(P_y))
end

MI(; context...) = hcubature(((x,y),)->observed_information(x,y; context...), [-10θ, -1-10σ], [10θ, 1+10σ]; reltol=ε)[1]

ε = 1e-9
σ = 0.1
θ = 1.0

f(x; slope=2.0, context...) = tanh(slope*x)
f⁻(ŷ; slope=2.0, context...) = atanh(ŷ)/slope
f⁻′(ŷ; slope=2.0, context...) = 1.0/(slope*(1.0-ŷ^2))

function test_mis()
    mi = zeros(100)
    slopes = 10 .^ LinRange(-2, 2, length(mi))
    for (i,slope) ∈ enumerate(slopes)
        mi[i] = MI(θ=θ, σ=σ, slope=slope)
        # println("Got MI = $(mi[i]) for slope=$(slope).")
    end
    return (slopes=slopes, mi=mi)
end

(slopes,mi) = test_mis()
slope_opt = slopes[argmax(mi)]

mi_fun = LinearInterpolation(slopes, mi)

push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\usetikzlibrary{backgrounds}")


xlims_plt = -3, 3
ylims_plt = -1.1, 1.1
x_plt = LinRange(xlims_plt...,200)
y_plt = LinRange(-0.995, 0.995, 200)
y_plt_prior = LinRange(-1.25, 1.25, 200)
y_true = [-0.90, 0.0, 0.90]
special_slopes = [10^-1, slope_opt, 10^1]

fig = begin
    fig = @pgf TikzPicture()
    
    ax = @pgf Axis(    
        {
            axis_background_style={fill="none"},
            width=raw"\textwidth-1cm",
            height=raw"4cm",
            xlabel = raw"$\textrm{gain }f'(0)$",
            ylabel = raw"$\mathcal{I}(X;Y)$",
            xmode = "log",
            clip="false",
            xmin="10^-2",
            xmax="10^2",
            no_markers,
            name="curve",
        },
        # draw curve
        PlotInc({line_width="2pt", smooth}, Table(slopes, mi)),
        # draw asymptotes
        Plot({color="black", line_width="1pt", dashed}, Table([slopes[1],slopes[end]], [0,0])),
        Plot({color="black", line_width="1pt", dashed}, Table([slopes[1],slopes[end]], [1,1])),            # draw vertical lines
        # # draw vertical lines
        # (VLine({ color="gray", line_width="1pt"}, slope) for slope ∈ special_slopes)...,
        ("\\node (c-0-$(i)) at ($(slope),$(mi_fun(slope))) {\$\\pgfuseplotmark{*}\$}; " for (i,slope) ∈ enumerate(special_slopes))...
    )

    push!(fig,ax)

    for (i,slope) ∈ enumerate(special_slopes)
        context = Dict(:slope=>slope, :σ=>σ, :θ=>θ)
        
        # draw function plot
        ax_f = @pgf Axis(
            {
                axis_background_style={fill="black!5!white", fill_opacity=0.5},
                height="4cm",
                width=raw"0.4*\textwidth",
                xtick=[-1,0,1], 
                ytick=[-1,0,1], 
                xlabel=raw"$X$", 
                ylabel=raw"$Y=f(X)+\eta$",
                xmin=xlims_plt[1],
                xmax=xlims_plt[2],
                ymin=ylims_plt[1],
                ymax=ylims_plt[2],
            },
            # draw function
            Plot(
                {
                    color="black",
                    line_width="2pt",
                },
                Table(x_plt, f.(x_plt; context...))
            ),
            # draw horizontal lines
            (HLine({ dashed, index_of_colormap=j, line_width="1pt", }, yy ) for (j,yy) ∈ enumerate(y_true))...,
            "\\coordinate (c-$(i)-1) at (axis description cs: 0,0);",
            "\\coordinate (c-$(i)-2) at (axis description cs: 1,0);",
            # add slope label
            "\\path (-1,-$(slope)) -- node[midway, sloped, anchor=south] {s=$(slope)} (1,$(slope));",
        )

        # draw x distributions
        ax_x = @pgf Axis(
            {
                at=[raw"($(curve.north west)+(0cm,3.5cm)$)",raw"($(curve.north)+(0cm,7cm)$)",raw"($(curve.north east)+(0cm,3.5cm)$)"][i],
                anchor=["south west", "south", "south east"][i],
                axis_background_style={fill="none"},
                height="2.5cm",
                width=raw"0.4*\textwidth",
                hide_axis,
                xlabel=raw"$p(X)^*, p(X|y_i)^*$", 
                xmin=xlims_plt[1],
                xmax=xlims_plt[2],
                ymin=0,
            },
            Plot(
                {
                    draw = "none",
                    fill="gray",
                    fill_opacity=0.5,
                },
                Table(x_plt,p_x.(x_plt; context...)./maximum(p_x.(x_plt; context...)))
            ),
            # draw posteriors
            (Plot(
                {
                    draw = "none",
                    index_of_colormap=j, 
                    fill=".",
                    fill_opacity=0.5
                },
                Table([xlims_plt[1];x_plt;xlims_plt[2]], [0;p_x_y.(x_plt,yⱼ; context...)./maximum(p_x_y.(x_plt,yⱼ; context...));0])
            ) for (j,yⱼ) ∈ enumerate(y_true))...,
        )
        
        # draw y distributions
        ax_y = @pgf Axis(
            {
                axis_background_style={fill="none"},
                height="4cm",
                width="2.5cm",
                hide_axis,
                xlabel=raw"$p(Y)^*,p(f(X)|y_i)^*$", 
                xmin=0,
                ymin=ylims_plt[1],
                ymax=ylims_plt[2],
            },
            # draw prior
            Plot(
                {
                    draw = "none",
                    fill="gray",
                    fill_opacity=0.5,
                },
                Table([0.0; p_y.(y_plt; context...)./maximum(p_y.(y_plt; context...)); 0.0], [-1.0; y_plt; 1.0])
            ),
            # # draw posteriors
            # (Plot(
            #     {
            #         draw = "none",
            #         index_of_colormap=j, 
            #         fill=".",
            #         fill_opacity=0.5
            #     },
            #     Table([0.0; p_ŷ_y.(y_plt,yⱼ; context...)./maximum(p_ŷ_y.(y_plt,yⱼ; context...)); 0], [-1.01; y_plt; 1.01])
            # ) for (j,yⱼ) ∈ enumerate(y_true))...,
        )

        gp = @pgf GroupPlot(
            {
                group_style = {
                    group_name="group$(i)",
                    group_size="2 by 2",
                    horizontal_sep="0pt",
                    vertical_sep="0pt",
                    xticklabels_at="edge bottom",
                    yticklabels_at="edge left"
                },
            },
            ax_x, {"group/empty plot",width="2.5cm",height="2.5cm", axis_background_style={fill="none"}}, ax_f, ax_y
        )
        #     plt[xidx][:yaxis][:flip]=true
        #     plt[fidx][:xaxis][:mirror]=true
        #     plt[yidx][:yaxis][:mirror]=true
        #     plt[yidx][:xaxis][:mirror]=true
        #     plt[eidx][:yaxis][:mirror]=true
        push!(fig, gp)
    end

    push!(fig,raw"\begin{scope}[on background layer]",("\\path[shade,bottom color=black!50!white,top color=black!10!white, fill opacity=0.2] (c-0-$(i).center) -- (c-$(i)-1) -- (c-$(i)-2) -- cycle;" for i ∈ 1:3)...,raw"\end{scope}")

    return fig
end
