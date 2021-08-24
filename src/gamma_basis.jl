using OrdinaryDiffEq
a=100.0
function df(du,u,p,t) 
	du[1:end-1] .= a*(u[2:end] .- u[1:end-1])
	du[end] = -a*u[end]
end

u0 = [0.0,0.0,0.0,0.0,1.0]
tspan = (0,0.1)
prob = ODEProblem(df, u0, tspan)
sol = solve(prob, Tsit5(), saveat=0.001)

leg=["\$1^\\text{st}\$ tap", "\$2^\\text{nd}\$ tap", "\$3^\\text{rd}\$ tap", "\$4^\\text{th}\$ tap", "\$5^\\text{th}\$ tap"]

plts = []
for i in 1:5
	push!(plts, @pgf PlotInc(
			{ no_marks, very_thick },
			Table(; x = sol.t, y = sol[end+1-i,:])
		)
	)
	push!(plts, @pgf LegendEntry(leg[i]))
end

p = @pgf TikzPicture(
	Axis({
		height="4cm", 
		width=raw"\textwidth",
		# "scale only axis", 
		xmin=0,
		xmax=0.1,
		ymin=0,
		ymax=1,
		xlabel="time in seconds",
		ylabel="impulse response",
		"axis background/.style"={fill="none"},
		},
		plts...
	)
)