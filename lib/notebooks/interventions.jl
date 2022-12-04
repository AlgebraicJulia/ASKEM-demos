### A Pluto.jl notebook ###
# v0.19.16

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 82d23846-7045-11ed-152b-49bacb843ba4
begin
	using Pkg
	Pkg.activate("..")
	using Revise, ASKEM
	using DifferentialEquations
	using Plots

	using Catlab.Theories
	using Catlab.Graphics
	using Catlab.CategoricalAlgebra
	using Catlab.Programs.RelationalPrograms

	using AlgebraicPetri
	using AlgebraicPetri.Epidemiology

	using ModelingToolkit
	using PlutoUI
	using HypertextLiteral
end

# ╔═╡ 72daa7df-e8c2-4913-b445-9fbbd8b885ec
md"""
Policies:
- Implement restrictions at fixed date
- Implement restrictions after cases are above a given number

Potential real world challenges (might not get to these):
- Have to give people time to prepare for lockdown
- Reported cases lag behind actual cases
- Effects of lockdown on transmission unknown

Problem:
Make a toolkit for evaluating different policies in different scenarios. This should include explorative evaluation and probabilistic evaluation.

TODO:
Figure appropriate scaling for all variables to correspond to actual situation.
"""

# ╔═╡ b003c833-8d56-4a7d-af3d-7007328368ef
sir_uwd = @relation (s, i, r) where (s, i, r) begin
    infection(s, i)
    recovery(i, r)
end;

# ╔═╡ ede2718c-489d-42df-a39f-fdec4646a9c1
sir_petri = apex(oapply_epi(sir_uwd));

# ╔═╡ 9b02c29f-012a-4b90-9b36-97871f9b0075
md"""
TODO: Should have a selection of different Petri nets, SIR, SEIR, etc.

For each Petri net, we should expose the variables of interest, i.e. case counts and hospitalization, as expressions in the ModelingToolkit variables.
"""

# ╔═╡ f62b4be9-9da4-41c9-a6e4-7da5653403e7
@named sir = ASKEM.PetriMTK.rate_equation(sir_petri)

# ╔═╡ e821e0fe-16fe-4320-8ea7-2d4dd91f3cb9
@bind intervention_params Sliders(
	"Intervention Parameters",
	["case_threshhold" => 0:5.0:500, 
	 "lockdown_effect" => 0:0.01:1.0,
	 "time_of_lockdown" => 0:0.01:1.0
	]
)

# ╔═╡ e542286b-bfa8-43d9-b150-6476d5580825
(case_threshold, lockdown_effect, time_of_lockdown) = intervention_params

# ╔═╡ 620e4c0e-b5af-412b-8a96-b758cb5387f9
@bind policy_choice Select([
	:case => "Caseload Triggered Lockdown",
	:time => "Time Triggered Lockdown"
])

# ╔═╡ ceda9aac-6e23-4606-9233-d76758d0913f
policy = Dict(
	:case => caseload_triggered_lockdown(sir, sir.I, case_threshold, lockdown_effect),
	:time => time_triggered_lockdown(sir, time_of_lockdown, lockdown_effect)
)[policy_choice];

# ╔═╡ b895c52e-1592-4905-a38c-7c9d70b92509
@named sir_intervened = ModelingToolkit.compose(policy, sir);

# ╔═╡ bbc73193-584e-4681-b30f-e6b19438eae7
@bind init Sliders("Initial Concentrations", ["S" => 0:50.0:1000, "I" => 0:5.0:100, "R" => 0:5.0:100])

# ╔═╡ c6c820af-1a8c-46c2-827e-547b9ab817c5
@bind rates Sliders("Rates", ["inf" => 0:0.001:0.1, "rec" => 0:0.1:2.0])

# ╔═╡ fce68171-13a9-4af3-a012-8d8b1ce763c6
prob = ODEProblem(sir_intervened, [init..., rates...], (0.0, 2.0), []);

# ╔═╡ 66f0ce87-2c08-4ad7-838b-56f9ed5029ac
plot(solve(prob))

# ╔═╡ 7e16d79d-1d1b-4003-9db2-de8de296a40f
md"""
Thoughts about the SIR model:

As long as I > 0 and inf > 0, the SIR model always says that S -> 0 eventually. Thus, the SIR model projects that lockdowns are ultimately useless except for lengthening the period over which people are infected. It does not account for the possibility of local eradication, where "burned over patches" of people who have been sick "firestop" transmission into some communities who never have covid introduced to them, which is a mechanism used to attempt to partially explain why covid infections come in waves.

Modeling this will need an agent-based, spatial model, and I expect that we might get interestingly different results from this.
"""

# ╔═╡ Cell order:
# ╠═82d23846-7045-11ed-152b-49bacb843ba4
# ╟─72daa7df-e8c2-4913-b445-9fbbd8b885ec
# ╠═b003c833-8d56-4a7d-af3d-7007328368ef
# ╠═ede2718c-489d-42df-a39f-fdec4646a9c1
# ╟─9b02c29f-012a-4b90-9b36-97871f9b0075
# ╠═f62b4be9-9da4-41c9-a6e4-7da5653403e7
# ╟─e821e0fe-16fe-4320-8ea7-2d4dd91f3cb9
# ╟─e542286b-bfa8-43d9-b150-6476d5580825
# ╟─620e4c0e-b5af-412b-8a96-b758cb5387f9
# ╠═ceda9aac-6e23-4606-9233-d76758d0913f
# ╠═b895c52e-1592-4905-a38c-7c9d70b92509
# ╟─bbc73193-584e-4681-b30f-e6b19438eae7
# ╟─c6c820af-1a8c-46c2-827e-547b9ab817c5
# ╠═fce68171-13a9-4af3-a012-8d8b1ce763c6
# ╠═66f0ce87-2c08-4ad7-838b-56f9ed5029ac
# ╟─7e16d79d-1d1b-4003-9db2-de8de296a40f
