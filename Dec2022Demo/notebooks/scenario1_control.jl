### A Pluto.jl notebook ###
# v0.19.12

using Markdown
using InteractiveUtils

# ╔═╡ 9aec393e-a083-45f5-ad73-e6bef22bb056
begin 
	using Pkg#, Revise
	Pkg.activate(Base.current_project())
	
	using ModelingToolkit
	using AlgebraicPetri
	using AlgebraicPetri: Graph
    using AlgebraicPetri.Epidemiology
	using AlgebraicPetri.BilayerNetworks

	using Catlab, Catlab.Theories
	using Catlab.CategoricalAlgebra
	using Catlab.Graphics
	using Catlab.Graphics: Graphviz
	import Catlab.CategoricalAlgebra: migrate!
	using Catlab.WiringDiagrams
	using Catlab.Programs
	using Catlab.Programs.RelationalPrograms
	import Catlab.WiringDiagrams.DirectedWiringDiagrams: WiringDiagramACSet
	import Catlab.CategoricalAlgebra.CSets: parse_json_acset

	using OrdinaryDiffEq
	using Optimization
	using OptimizationOptimisers

	import ASKEM.Oct2022Demo: sumvarsbyname, MakeReactionSystem
	using ASKEM.Upstream: presentationToLabelledPetriNet, deserialize_wiringdiagram, vectorfield
	using ASKEM.Dec2022Demo: formSIRD, formTVParams, solveODE, zeroVal, runControlOptim, makeK, runControlAuto, draw, sig, invsig

	using Plots
end

# ╔═╡ 1df5e4a8-e761-4767-aed4-866534398922
md"""## Scenario 1: Keeping Hospitalization Below Threshold"""

# ╔═╡ 98ba5842-8ff9-4bce-a0a1-04c228b9ff9c
md"""We have limited data, so we interpret the scenario as:
- SIRD Model
- H is a function (fixed fraction) of I
- using a set of parameters that best generate example dynamics and demonstrate control
"""

# ╔═╡ 943d7c37-ade5-4a4d-ae0b-e2bf390fcf6f
md"""### SIRD Model (with H Observed)"""

# ╔═╡ 85ddcf05-c9c5-4436-8511-3713dbbe7694
begin
	SIRD = read_json_acset(LabelledPetriNet,"../SIRD.json")
	AlgebraicPetri.Graph(SIRD)
end

# ╔═╡ 59b234fe-0376-497c-bea3-997ab50422c5
md"""H = hosp_rt * I"""

# ╔═╡ c7535760-7c2d-4552-a39a-7b51af04a92f
md"""### Set Up Example"""

# ╔═╡ cf98d2fa-2946-46ff-ba5a-00f8ee9d263c
begin
	u0 = [999,1,0,0]
	p_fixed = [0.000025,0.005,0.001]
	tspan = (1,1000)
	alpha_policy = invsig(0.05)
	tstart_policy = 1
	hosp_rt = 0.5
	thresh_H = 125
	t_start = 200
	alpha_init = [0.0]
end

# ╔═╡ 61a6b383-b9e2-4acf-8bb7-0fa926a0ec12
md"""### Control Using Fixed Policy"""

# ╔═╡ df1e0fa3-dbea-4a61-a593-9ef6c517f90b
md"""- Run with no intervention
- Run with 5% decrease starting beginning of time span"""

# ╔═╡ b323ac71-837a-4da1-ac81-05ac1ca9d600
begin
	ControlWorkflow = read_json_acset(LabelledPetriNet,"../s1_cntrl_wf_present.json")
	AlgebraicPetri.Graph(ControlWorkflow)
end

# ╔═╡ 6762c0ea-a453-451c-80c1-566a1389a45c
begin
	s1_sird_cntrl_policy = deserialize_wiringdiagram("../s1_sird_cntrl_policy.json")
	draw(s1_sird_cntrl_policy)
end

# ╔═╡ be745863-ca88-49ca-bc62-553cbafaefb7
md"""#### Without Intervention"""

# ╔═╡ d797cea2-7c16-499d-8659-4ff352f8c029
formSIRD()

# ╔═╡ 4b2a1018-c683-401e-8da2-e6792c5c162b
begin
	free_hom_expr = to_hom_expr(FreeBiproductCategory,s1_sird_cntrl_policy)
	free_jfunc = Catlab.Programs.GenerateJuliaPrograms.compile_expr(free_hom_expr)
end

# ╔═╡ 3fed27e7-e09c-4099-a04f-62ebf3ddac03


# ╔═╡ 7ae8fdc8-1111-46ea-bfcc-85cbc11e73c3
free_tv_sol = eval(free_jfunc)(u0,p_fixed,tspan,0,tstart_policy)

# ╔═╡ 8fddc83d-8993-4b82-b2e2-d833c8bd9f26
md"""#### Fixed policy - 5% decrease in infection rate on first day"""

# ╔═╡ 2f818299-058c-4e0a-ad91-09cf4899c70e
begin 
	policy_hom_expr = to_hom_expr(FreeBiproductCategory,s1_sird_cntrl_policy)
	policy_jfunc = Catlab.Programs.GenerateJuliaPrograms.compile(policy_hom_expr)
end

# ╔═╡ 28d7e1e3-4ecd-4c8c-9ddc-7cb6510c28b8
policy_tv_sol = policy_jfunc(u0,p_fixed,tspan,alpha_policy,tstart_policy)

# ╔═╡ 505cd6b7-f26c-42e0-812b-6962255d3648
md"""### Control via optimization - keep hospitalizations under 125 w/ decrease starting day 200"""

# ╔═╡ f5a2b608-d8ea-4cfa-bc53-80ee170443eb
begin
	s1_sird_cntrl_optim = deserialize_wiringdiagram("../s1_sird_cntrl_optim.json")
	draw(s1_sird_cntrl_optim)
end

# ╔═╡ 4d731c25-b931-4f70-96e7-bbe231b29db0
begin
	optim_hom_expr = to_hom_expr(FreeBiproductCategory,s1_sird_cntrl_optim)
	optim_jfunc = Catlab.Programs.GenerateJuliaPrograms.compile(optim_hom_expr)
end

# ╔═╡ 2f16a5a7-8ca5-43bc-aa29-1436d4e68fe0
alpha, optim_tv_sol, t, obs_hosp = 	optim_jfunc(u0,p_fixed,tspan,hosp_rt,thresh_H,t_start,alpha_init)

# ╔═╡ 7366703c-7976-4cc0-ab6e-da39054270d9
md"""### Automated LQR control"""

# ╔═╡ bb48c3e5-5ba0-47f0-8df7-a4a90863dad9
begin
	s1_sird_cntrl_auto = deserialize_wiringdiagram("../s1_sird_cntrl_auto.json")
	draw(s1_sird_cntrl_auto)
end

# ╔═╡ dc8b5de5-721b-4f7f-ad2d-c9e0a077d18e
begin
	auto_hom_expr = to_hom_expr(FreeBiproductCategory,s1_sird_cntrl_auto)
	auto_jfunc = Catlab.Programs.GenerateJuliaPrograms.compile(auto_hom_expr)
end

# ╔═╡ 7f5c4e33-1c8b-41a1-ba73-52d3643fa008
auto_tv_sol, t_auto = auto_jfunc(u0,p_fixed,tspan)

# ╔═╡ Cell order:
# ╠═9aec393e-a083-45f5-ad73-e6bef22bb056
# ╟─1df5e4a8-e761-4767-aed4-866534398922
# ╠═98ba5842-8ff9-4bce-a0a1-04c228b9ff9c
# ╠═943d7c37-ade5-4a4d-ae0b-e2bf390fcf6f
# ╠═85ddcf05-c9c5-4436-8511-3713dbbe7694
# ╠═59b234fe-0376-497c-bea3-997ab50422c5
# ╠═c7535760-7c2d-4552-a39a-7b51af04a92f
# ╠═cf98d2fa-2946-46ff-ba5a-00f8ee9d263c
# ╠═61a6b383-b9e2-4acf-8bb7-0fa926a0ec12
# ╠═df1e0fa3-dbea-4a61-a593-9ef6c517f90b
# ╠═b323ac71-837a-4da1-ac81-05ac1ca9d600
# ╠═6762c0ea-a453-451c-80c1-566a1389a45c
# ╠═be745863-ca88-49ca-bc62-553cbafaefb7
# ╠═d797cea2-7c16-499d-8659-4ff352f8c029
# ╠═4b2a1018-c683-401e-8da2-e6792c5c162b
# ╠═3fed27e7-e09c-4099-a04f-62ebf3ddac03
# ╠═7ae8fdc8-1111-46ea-bfcc-85cbc11e73c3
# ╠═8fddc83d-8993-4b82-b2e2-d833c8bd9f26
# ╠═2f818299-058c-4e0a-ad91-09cf4899c70e
# ╠═28d7e1e3-4ecd-4c8c-9ddc-7cb6510c28b8
# ╠═505cd6b7-f26c-42e0-812b-6962255d3648
# ╠═f5a2b608-d8ea-4cfa-bc53-80ee170443eb
# ╠═4d731c25-b931-4f70-96e7-bbe231b29db0
# ╠═2f16a5a7-8ca5-43bc-aa29-1436d4e68fe0
# ╠═7366703c-7976-4cc0-ab6e-da39054270d9
# ╠═bb48c3e5-5ba0-47f0-8df7-a4a90863dad9
# ╠═dc8b5de5-721b-4f7f-ad2d-c9e0a077d18e
# ╠═7f5c4e33-1c8b-41a1-ba73-52d3643fa008
