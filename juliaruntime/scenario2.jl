### A Pluto.jl notebook ###
# v0.19.16

using Markdown
using InteractiveUtils

# ╔═╡ bd4b4c28-2fc0-4dcf-93d6-32af8310044d
# ╠═╡ show_logs = false
let
	import Pkg
	Pkg.activate(".")
    include("lib.jl")
	nothing
end

# ╔═╡ d0974f57-02a0-421b-9f2b-5088a7ccc01d
using CSV, DataFrames, Plots, StatsPlots

# ╔═╡ 4d8041ea-a483-4cac-87c9-b85cd64f5ab2
using AxisArrays

# ╔═╡ c9649f10-bc07-49d5-b96e-41271216dbae
md"""### Obtain sim plan and parameters from front end"""

# ╔═╡ e23fd17c-8c7c-4090-9a27-772e8260b615
run(`cat forecast_plan.json`);

# ╔═╡ ec6be9d3-bad2-47ce-ab50-00d83b8245a4
state, params, tspan = [0.9,0.5,0,0], [0.2, 0.5, 0.01], (0,20)

# ╔═╡ d771ae6d-632d-4332-a8f9-0ba63495bb56
md"""### Execute and return result"""

# ╔═╡ 1bd8893e-120d-41e4-aeaa-a4d88ea323f5
forecast = deserialize_wiringdiagram("forecast_plan.json")

# ╔═╡ 99926aa8-2de8-4884-aeae-bd8351a6d603
begin
	sim = interpreter(forecast, (; formSIRD, mtk_simulate, soln_to_csv))
	sim.set_inputs!(state, params, tspan)
	out = sim.get_outputs!()[1]
	sim.stop!()
	println(out)
end

# ╔═╡ 6a96b960-22a4-4395-9549-d29a01844af1
df = CSV.read(IOBuffer(out), DataFrame)

# ╔═╡ e8f91749-19a9-4722-8426-66ccc222b46b
@eval @df df plot(:timestamp, [:var"S(t)", :var"I(t)", :var"R(t)", :var"D(t)"], label=["S" "I" "R" "D"])

# ╔═╡ e62b9b86-708f-441d-9574-c71efbb6d614
md"""### Get stratified model from front end"""

# ╔═╡ bebf62d7-0d1a-489f-8bd3-b22cb77c86ee
stratified_mdl_json = read("../Dec2022Demo/outputs/sird_age7_vax.json",String);

# ╔═╡ 3f4570ca-74d0-451a-bc4f-51d82f051a47
println(stratified_mdl_json)

# ╔═╡ 5458f1db-7dad-4dae-b8a0-11fda964aa0f
run(`cat forecast_plan2.json`)

# ╔═╡ 4b3a2297-20f9-4c4e-8dfd-c2d1fdee9046
state2, params2 = rand(56), rand(217)

# ╔═╡ 8c9b2eaf-e802-4444-8462-4065fa632147
md"""### Again, execute"""

# ╔═╡ a8c3f724-e14b-42cd-9c1f-42e28f63bff4
forecast2 = deserialize_wiringdiagram("forecast_plan2.json")

# ╔═╡ 67b56e7b-bd72-46c3-b138-5d72a354c7fc
begin
	sim2 = interpreter(forecast2, (; json_to_mdl, mtk_simulate, soln_to_csv))
	sim2.set_inputs!(stratified_mdl_json, state2, params2, tspan)
	out2 = sim2.get_outputs!()[1]
	sim2.stop!()
	println(out2)
end

# ╔═╡ 0b783a21-c636-42eb-8f58-38381a5531a3
_df2 = CSV.read(IOBuffer(out2), DataFrame);

# ╔═╡ 0bccdaf1-3666-492a-9c86-b36ff05fc8ed
begin
	function mapname(sym)
	    s = string(sym)
	    Symbol('"' in s ? join(split(s, '"')[2:2:end], '-') : s)
    end
    df2 = DataFrame(mapname.(names(_df2)) .=> eachcol(_df2))
end

# ╔═╡ 4d041d1d-6dfc-430b-8a76-8f29cfa29c47
length(eachcol(df2))

# ╔═╡ fdd5dc1f-9bd2-4dfa-9848-6f9b7d7df59a
md"""### Form data cube for querying"""

# ╔═╡ 43fc104d-431f-4e6f-9fb7-02ba798af98a
cube = AxisArray(zeros(length(eachrow(df2)), 7, 4, 2),
	       time = df2[:, :timestamp],
	       age = [ Symbol("Age$i") for i=1:7],
	       state = [:S, :I, :R, :D],
	       v = [:U, :V]);

# ╔═╡ 24f6f0bb-51cd-4fa1-8ec2-9f7c25fe1daa
for c = 2:57
	s, a, v = Symbol.(split(names(df2)[c], '-'))
	cube[:, a, s, v] = df2[:, c]
end

# ╔═╡ 9a3546e7-70d7-4945-990d-198c3b3325b1
cube

# ╔═╡ 20e5dc92-44a9-487c-b6e4-0e8001d5e4af
cube[0.0..0.3, :Age2, :R, :U]

# ╔═╡ 0c2500d1-647d-46e7-981c-13f08c65ec83
# TODO: output as netcdf for querying in other tools

# ╔═╡ Cell order:
# ╠═bd4b4c28-2fc0-4dcf-93d6-32af8310044d
# ╟─c9649f10-bc07-49d5-b96e-41271216dbae
# ╠═e23fd17c-8c7c-4090-9a27-772e8260b615
# ╠═ec6be9d3-bad2-47ce-ab50-00d83b8245a4
# ╟─d771ae6d-632d-4332-a8f9-0ba63495bb56
# ╠═1bd8893e-120d-41e4-aeaa-a4d88ea323f5
# ╠═99926aa8-2de8-4884-aeae-bd8351a6d603
# ╠═d0974f57-02a0-421b-9f2b-5088a7ccc01d
# ╠═6a96b960-22a4-4395-9549-d29a01844af1
# ╠═e8f91749-19a9-4722-8426-66ccc222b46b
# ╟─e62b9b86-708f-441d-9574-c71efbb6d614
# ╠═bebf62d7-0d1a-489f-8bd3-b22cb77c86ee
# ╠═3f4570ca-74d0-451a-bc4f-51d82f051a47
# ╠═5458f1db-7dad-4dae-b8a0-11fda964aa0f
# ╠═4b3a2297-20f9-4c4e-8dfd-c2d1fdee9046
# ╟─8c9b2eaf-e802-4444-8462-4065fa632147
# ╠═a8c3f724-e14b-42cd-9c1f-42e28f63bff4
# ╠═67b56e7b-bd72-46c3-b138-5d72a354c7fc
# ╠═0b783a21-c636-42eb-8f58-38381a5531a3
# ╠═0bccdaf1-3666-492a-9c86-b36ff05fc8ed
# ╠═4d041d1d-6dfc-430b-8a76-8f29cfa29c47
# ╠═4d8041ea-a483-4cac-87c9-b85cd64f5ab2
# ╟─fdd5dc1f-9bd2-4dfa-9848-6f9b7d7df59a
# ╠═43fc104d-431f-4e6f-9fb7-02ba798af98a
# ╠═24f6f0bb-51cd-4fa1-8ec2-9f7c25fe1daa
# ╠═9a3546e7-70d7-4945-990d-198c3b3325b1
# ╠═20e5dc92-44a9-487c-b6e4-0e8001d5e4af
# ╠═0c2500d1-647d-46e7-981c-13f08c65ec83
