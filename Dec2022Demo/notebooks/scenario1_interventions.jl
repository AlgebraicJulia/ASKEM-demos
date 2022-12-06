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
This notebook allows you to choose between several epidemiological Petri nets and several different policies, and then simulate what the effect of those policies would be.
"""

# ╔═╡ 82c9f255-41c2-4933-b31b-6558f9243495
md"""
First, select a Petri net.
"""

# ╔═╡ 62b804b2-d355-42fb-8e4c-f5eb9eed1696
petri_options = [
	EpidemiologicalPetriNets.SIR => "SIR",
	EpidemiologicalPetriNets.SIRD => "SIRD",
	EpidemiologicalPetriNets.SEIR => "SEIR",
	EpidemiologicalPetriNets.SEIRD => "SEIRD"
];

# ╔═╡ 9ea0277f-ed20-44e8-bbf2-16e8815ee65b
@bind petri Select(petri_options)

# ╔═╡ f62b4be9-9da4-41c9-a6e4-7da5653403e7
@named model = ASKEM.PetriMTK.rate_equation(petri)

# ╔═╡ a437c4ce-f5ae-43bb-ac44-c14141fe368a
md"""
Then choose a policy. The two choices are a lockdown at a fixed, predefined time, and a lockdown triggered by a certain number of cases.
"""

# ╔═╡ 620e4c0e-b5af-412b-8a96-b758cb5387f9
@bind policy_choice Select([
	:case => "Caseload Triggered Lockdown",
	:time => "Time Triggered Lockdown"
])

# ╔═╡ f1f9107d-836f-4807-98d0-08a96e9b6a5e
md"""
The parameters for the policy can then be set with the following sliders.
"""

# ╔═╡ e821e0fe-16fe-4320-8ea7-2d4dd91f3cb9
@bind intervention_params Sliders(
	"Intervention Parameters",
	["case_threshhold" => 0:0.01:1.0, 
	 "lockdown_effect" => 0:0.01:1.0,
	 "time_of_lockdown" => 0:1.0:50.0
	]
)

# ╔═╡ e542286b-bfa8-43d9-b150-6476d5580825
(case_threshold, lockdown_effect, time_of_lockdown) = intervention_params

# ╔═╡ ceda9aac-6e23-4606-9233-d76758d0913f
policy = Dict(
	:case => caseload_triggered_lockdown(model, model.I, case_threshold, lockdown_effect),
	:time => time_triggered_lockdown(model, time_of_lockdown, lockdown_effect)
)[policy_choice];

# ╔═╡ b895c52e-1592-4905-a38c-7c9d70b92509
@named model_intervened = ModelingToolkit.compose(policy, model);

# ╔═╡ bbc73193-584e-4681-b30f-e6b19438eae7
@bind init Sliders("Initial Concentrations", [string(sname(petri, i)) => 0:0.001:1 for i in 1:ns(petri)])

# ╔═╡ c6c820af-1a8c-46c2-827e-547b9ab817c5
@bind rates Sliders("Rates", [string(tname(petri, i)) => 0:0.001:1 for i in 1:nt(petri)])

# ╔═╡ fce68171-13a9-4af3-a012-8d8b1ce763c6
prob = ODEProblem(model_intervened, [init..., rates...], (0.0, 50.0), []);

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
# ╟─82c9f255-41c2-4933-b31b-6558f9243495
# ╟─62b804b2-d355-42fb-8e4c-f5eb9eed1696
# ╠═9ea0277f-ed20-44e8-bbf2-16e8815ee65b
# ╠═f62b4be9-9da4-41c9-a6e4-7da5653403e7
# ╟─a437c4ce-f5ae-43bb-ac44-c14141fe368a
# ╟─620e4c0e-b5af-412b-8a96-b758cb5387f9
# ╟─f1f9107d-836f-4807-98d0-08a96e9b6a5e
# ╠═e821e0fe-16fe-4320-8ea7-2d4dd91f3cb9
# ╟─e542286b-bfa8-43d9-b150-6476d5580825
# ╟─ceda9aac-6e23-4606-9233-d76758d0913f
# ╠═b895c52e-1592-4905-a38c-7c9d70b92509
# ╠═bbc73193-584e-4681-b30f-e6b19438eae7
# ╠═c6c820af-1a8c-46c2-827e-547b9ab817c5
# ╠═fce68171-13a9-4af3-a012-8d8b1ce763c6
# ╠═66f0ce87-2c08-4ad7-838b-56f9ed5029ac
# ╟─7e16d79d-1d1b-4003-9db2-de8de296a40f
