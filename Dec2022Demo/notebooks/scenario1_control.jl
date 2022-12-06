### A Pluto.jl notebook ###
# v0.19.12

using Markdown
using InteractiveUtils

# ╔═╡ 915f4df0-bc3a-11ec-2560-c9772e143679
begin 
	using Pkg#, Revise
	Pkg.activate(Base.current_project())
	# Pkg.instantiate()
	using Catlab.CategoricalAlgebra
	using Catlab.Present, Catlab.Theories
	using AlgebraicPetri
	using AlgebraicPetri: Graph
    # using ModelExploration.Explore
	using Plots
	# using Test

end;

# ╔═╡ aad0407a-4979-4ea4-94d3-38f1bb6d3de1
begin 
	using Random
	# rng = Random.default_rng()
	Random.seed!(1234)
end;

# ╔═╡ 60d7decb-a7fa-494e-93ca-26d8b957cfcf
include("LibForOct2022Demo.jl");

# ╔═╡ 7603c764-9e3f-4a66-9ee6-e7694e5070a3
include("LibForDec2022Demo.jl")

# ╔═╡ 9175d027-b5eb-4676-bed3-e084a882e102
# include("StratPetris.jl")

# ╔═╡ 943d7c37-ade5-4a4d-ae0b-e2bf390fcf6f
md"""## Control of SIR Model (with H Observed)"""

# ╔═╡ 85ddcf05-c9c5-4436-8511-3713dbbe7694
SIRD = LabelledPetriNet([:S, :I, :R, :D],
    :inf => ((:S,:I) => (:I,:I)),
    :recover => (:I=>:R),
    :death => (:I=>:D));

# ╔═╡ 59b234fe-0376-497c-bea3-997ab50422c5
AlgebraicPetri.Graph(SIRD)

# ╔═╡ 61a6b383-b9e2-4acf-8bb7-0fa926a0ec12
md"""## Set up ODEProblem with time-varying parameter, lagged h observation"""

# ╔═╡ 55e42e79-b7f1-4793-aaa7-38ecfe0f3b66
begin
	p_i(t) = 0.1+0.01*cos(t/100*4)
	p_r(t) = 0.05+0.01*cos(t/100*4+pi/4)
	p_d(t) = 0.01+0.002*cos(t/100*4+pi/2)
	p_t = [p_i,p_r,p_d]
end

# ╔═╡ 88313435-2633-40e3-94b1-661c15be2b30
begin
	tv_rxn = vectorfield(SIRD)
	tspan = (0,100)
	u0 = [999,1,0,0]
end

# ╔═╡ a984a49c-8fca-435c-a49b-09e72d5f6c48
tv_prob = ODEProblem(tv_rxn, u0, tspan, p_t)

# ╔═╡ c7a089e6-9cb1-4a18-a65d-ff4044204eec
tv_sol = solve(tv_prob, Tsit5())

# ╔═╡ ba40741a-3c3f-4e22-9ad1-4c150db5deb8
p_i

# ╔═╡ 8058af4a-4df7-45e3-82ff-53bf10908baf
tnames(SIRD)

# ╔═╡ 1e2aeb26-5f62-4f06-956a-1926befb923a
test = [0.1+0.01*cos(t/100*4),0.05+0.01*cos(t/100*4+pi/4),0.01+0.002*cos(t/100*4+pi/2)]

# ╔═╡ b9c09401-186e-43c3-98db-f36538a88616
test(2)

# ╔═╡ bcd3233e-2e0a-44f0-8f18-9123a9a89012
md"""## Simulate with steps"""

# ╔═╡ 323a8caf-bb26-4c7b-b6b3-524351a80e5c
md"""## Optimize min-rate to keep below"""

# ╔═╡ 6e7cca77-de89-42c5-a8ce-c016eb327fe3
md"""## Controller version, with two lags - change in rate from controller and lag in H"""

# ╔═╡ cfc6c757-8893-4d9c-add8-4fd6906ae156
md"""## Model Calibration of MIRA-Two-City Model"""

# ╔═╡ 09444324-aa7b-4e63-af94-fa8950ece8b6
md"""Here we simulate data from the stratified model (formed from the MIRA disease model and Semagrams two-city stratification model) and fit the model to the data sample."""

# ╔═╡ 287a44bb-e4d9-4a8b-b383-13d9c01e6e1e
md"""### Extract MIRA Model and Stratification"""

# ╔═╡ f0f007c1-6528-4a85-a007-b31881297f6a
md"""**Load Disease Model from Curated MIRA JSON**"""

# ╔═╡ c7a6ed94-e97a-4a96-8f07-1370e71038b9
begin
	@present TheoryMIRANet <: TheoryLabelledReactionNet begin
    	MID::AttrType
    	MCTX::AttrType
    	Template::AttrType
    	mira_ids::Attr(S, MID)
    	mira_context::Attr(S, MCTX)
    	template_type::Attr(T, Template)
    	parameter_name::Attr(T, Name)
    	# parameter_value::Attr(T, Rate)
	end
	@abstract_acset_type AbstractMIRANet <: AbstractLabelledReactionNet
	@acset_type MIRANet(TheoryMIRANet) <: AbstractMIRANet
end;

# ╔═╡ e2c11399-7878-4602-a276-e190857b3fa6
mdl_disease_mira = load_mira_curated("BIOMD0000000971_petri_curated.json", 0.05);

# ╔═╡ 6dec59bf-1fa0-4ad6-afdd-26ec2b92230a
begin
	filled_rates = mdl_disease_mira[:rate]
	mdl_disease_rxnet = LabelledReactionNet{Any,Any}()
	copy_parts!(mdl_disease_rxnet,mdl_disease_mira)
end;

# ╔═╡ 560fda7a-4f4e-49e9-bf09-d82557286a83
mdl_disease = LabelledPetriNet(mdl_disease_rxnet);

# ╔═╡ 2c24347b-9c55-4e14-876a-ea8e2867eaa2
md"""**Define Type System and Typed Disease Model**"""

# ╔═╡ f0767520-04f6-4e97-a889-cf5e45d70b4c
begin
	types′ = LabelledPetriNet([:Pop],
    	:infect=>((:Pop, :Pop)=>(:Pop, :Pop)),
    	:disease=>(:Pop=>:Pop),
    	:strata=>(:Pop=>:Pop))
	types = map(types′, Name=name->nothing)
end;

# ╔═╡ 8493b4cb-12fa-4a1b-9125-ea59523a889b
begin 
	SIRD_typed = homomorphism(SIRD, types;
    	initial=(T=[1,2,2],I=[1,2,3,3],O=[1,2,3,3]),
    	type_components=(Name=x->nothing,))
	@assert is_natural(SIRD_typed)
end

# ╔═╡ a3eca5f4-21c0-4147-b7d0-19ba3fd7a776
begin 
	mdl_disease_typed = homomorphism(mdl_disease, types;
    	initial=(T=[1,1,2,2,2,2,2,2,2,2,2,1,2,1,1,1],),
    	type_components=(Name=x->nothing,))
	@assert is_natural(mdl_disease_typed)
end

# ╔═╡ 12a51952-34ab-47a3-a318-c4f58e3f3c12
md"""**Extract Two-City Stratification Model**"""

# ╔═╡ 2cc7237d-3a8b-4b86-bb08-6cba2e9e9531
begin
	mdl_strat = LabelledPetriNet([:City1,:City2],
	    :travel12 => ((:City1)=>(:City2)),
    	:travel21 => ((:City2)=>(:City1)))
end;

# ╔═╡ 5b66805c-42ef-4611-a31b-962ba3f549ed
begin
	mdl_strat_typed = homomorphism(mdl_strat, types;
    	initial=(T=[3,3],), type_components=(Name=x->nothing,))
	@assert is_natural(mdl_strat_typed)
end;

# ╔═╡ 4653e385-beb8-49d6-a562-aae001db013e
 #=
 read_json_acset(StratPetriUntyped{Float64, Float64, String, String,Any}, "semagrams_stratified_model.json")
 read_json_acset(StratPetri,"semagrams_stratification_model.json")
 tc_ss = tostrataspec(tc_strat_petri)
 =#

# ╔═╡ c6720b8a-e709-4c47-944f-dd4e7d750cbb
md"""**Form Stratification**"""

# ╔═╡ 454f4cb1-1535-4aba-a824-181489d324de
begin
	disease_ss = StrataSpec(mdl_disease_typed, 
		[[:strata],[:strata],[:strata],[:strata],[],[],[:strata],[],[]])
	strat_ss = StrataSpec(mdl_strat_typed, [[:disease,:infect], 			[:disease,:infect]])
end;

# ╔═╡ 2f2d6daf-c0b8-4501-86a5-65a844c9359e
mdl_Mira_TC, obs_Mira_TC = stratify(disease_ss, strat_ss, types′);

# ╔═╡ e900b657-4898-4811-aa6c-e178d693e6bf
md"""For clarity we rename the disease and stratified models."""

# ╔═╡ 5aebef7b-6041-48b0-90cc-155975667b34
begin
	Mira = mdl_disease
	Mira_typed = mdl_disease_typed
	# mdl_Mira_TC = mdl_stratified 
	# obs_Mira_TC = obs_stratified
	TwoCity = mdl_strat
end;

# ╔═╡ bc5d6541-af3e-4101-a18e-fdcb8f3a7871
md"""For reference the MIRA disease model is"""

# ╔═╡ 1408b638-47ec-4cb7-9cd2-a2fea531eb8f
AlgebraicPetri.Graph(Mira)

# ╔═╡ 2f34272e-68d1-4356-ac7e-13aa4e221b7d
md"""And the MIRA-Two-City stratified model is"""

# ╔═╡ c04f6eaa-3c8a-4c8e-9311-63eb2e259315
AlgebraicPetri.Graph(mdl_Mira_TC)

# ╔═╡ 1d940676-91a4-4c58-a390-9fc276d706bc
md"""### Generate Data"""

# ╔═╡ d711e132-192c-47c1-b8a0-7070aa5f55f4
md"""**Define Observation Function**"""

# ╔═╡ a6e855e1-9992-4969-b56e-872866801a9b
md"""Crucially, to both fit the model to data and generate sample data from the model, we need to specify an **observation function** for the model, i.e., a function mapping the states to the observation space."""

# ╔═╡ be1f1295-a9d5-4dc8-9b6b-8269bbc15827
md"""This is necessary because, in general, not all of the states will be directly observable and/or the states may not be directly observed at all, being instead functions of some of the states."""

# ╔═╡ 03bf51d8-9735-458e-9d3a-816c9fcdde0c
md"""In particular, to perform model comparison and selection from fits to a given dataset, each model requires its own observation function into the given observation space.""" 

# ╔═╡ 91ad223d-c245-4aca-89ab-cd6ee8c6652d
md"""Resulting in the following stratified model"""

# ╔═╡ d0b5715e-1964-4957-884b-5315f2a61777
md"""Below, we explored two possible such observations functions."""

# ╔═╡ cb71ea78-9dc1-4bf9-871a-1a7f503f9371
md"""The first is the observation defined by the projection onto the states of the disease model, summing over the corresponding stratified states. This function is currently returned from the stratification process."""

# ╔═╡ 7b588922-125b-4866-868d-3ba9e63a307b
md"""The second is a restricted observation, selecting only a subset of states, specfically the I, H, and D states, as given by"""

# ╔═╡ 26f3d3ee-91ab-42f1-bc94-c8a7c689e1d2
function obs_IHD_from_func(obs_func, sol, sample_times)
	samples, labels = obs_func(sol, sample_times)
	indx(labels, l) = findall(isequal(l), labels)
	sumsamples(l) = sum(samples[:, indx(labels, l)], dims=2)
	i = sumsamples(:I)
	h = sumsamples(:H)
	d = sumsamples(:D)
	labels = [:I :H :D]
	# labels = reshape(["I", "H", "D"],1,3)

    return hcat(i, h, d), labels
end;

# ╔═╡ e2325e9e-94bc-490f-877b-0922757fbd62
md"""**Specify the initial states and rates**"""

# ╔═╡ 56bbf4fd-c7ff-4532-bfd3-12f802a9bbd5
md"""To generate sample data from the model, we provide rate paramters, initial states, and a time span and then simulate."""

# ╔═╡ 923fc804-c188-417e-9465-2380b139469b
md"""We start with 999 susceptible individuals and one infected in each city."""

# ╔═╡ 45811f40-ed90-4155-938a-f6fd91f7030c
begin
	true_mdl = mdl_Mira_TC
	# true_obs = (sol, times) -> obs_IHD_from_func(obs_Mira_TC, sol, times)
	true_obs = obs_Mira_TC
	# u0 = repeat([999.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0],2)
	u0_ = [999.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,
			1000.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
end;

# ╔═╡ da4e2714-7aec-4838-81a9-621c403130fb
md"""We use the rates of the component models to form the rates of the stratified model."""

# ╔═╡ ec05e159-50de-439a-8063-a9add8db5599
begin
	Mira_rates = filled_rates
	TC_rates = [0.2,0.2]
	function transidx(model, name)
		filter(parts(model, :T)) do i
			model[i, :tname] == name
		end
	end
	true_p = zeros(nt(true_mdl))
	for ii in 1:nt(true_mdl)
		if true_mdl[ii,:tname][1]==:strata
			true_p[ii] = TC_rates[transidx(TwoCity,true_mdl[ii,:tname][2])[1]] 
		else 
			true_p[ii] = Mira_rates[transidx(Mira,true_mdl[ii,:tname][1])[1]]
		end
	end
end;

# ╔═╡ 37a5cdeb-f98d-4eb0-9e5b-7e0935a26316
md"""**Generate Data Sample**"""

# ╔═╡ 1d3f03de-5b5b-4ea9-839c-a6ae67f72a22
begin
	tspan_ = (0.0,250.0)
	sample_data, sample_times, prob_real, true_sol, noiseless_data, data_labels = 		generateData(true_mdl, true_p, u0, tspan, 50, true_obs)
end;

# ╔═╡ 35a87a6e-4c70-400f-9e94-7b48e90e536f
plot(sample_data, labels=reshape(String.(data_labels),1,length(data_labels)))

# ╔═╡ 7e35a352-bc77-4eb1-9e5b-1a87fba90cc1
md"""### Calibrate Data"""

# ╔═╡ efa579a9-f8ee-4158-9ed5-1c17295a47c9
md"""To calibrate a model to sample data, we provide initial estimates for the rate paramters (in addition to the model and observation function)."""

# ╔═╡ a605d131-3e89-4701-a660-f32ad6182811
md"""Similar to the need for specifying an observation function in the generation of data, in general, model fitting requires specifying a loss function, e.g., weighting the loss from different observations differently."""

# ╔═╡ 39e2a785-e432-40a8-b679-49e9440cbf99
md"""Here, our calibration function offers the possibility to subselect states for inclusion in the loss evaluation. Currently, we are including all states from the disease model."""

# ╔═╡ 8956e904-e17f-4f9c-876f-8783dbd48bcd
begin
	# states_to_count = [:I,:H,:D]
	states_to_count = [:S,:E,:I,:A,:SQ,:H,:R,:EQ,:D]
	p_init = repeat([1.0e-1], nt(true_mdl))
	# p_init = true_p
	p_est, sol_est, loss = calibrate(true_mdl, true_obs, states_to_count, u0, p_init, 	sample_data, sample_times, data_labels)
end;

# ╔═╡ fc3b4271-4ea8-4e1f-a62f-69186d4fbac8
md"""### Results"""

# ╔═╡ 50ba48db-8466-4299-967e-825ab99fba9e
md"""**Estimated Parameters**"""

# ╔═╡ dcaf7bd1-c692-4067-a57c-d32f244a3054
hcat(true_mdl[:tname],p_est)

# ╔═╡ d65328b2-34f1-4321-b7d2-0dabca98e2e2
md"""**Estimated Observations**"""

# ╔═╡ 4d0e2692-cd10-455b-9030-7b59046504a9
plot_obs_w_ests(sample_times, sample_data, sol_est, true_obs)

# ╔═╡ 7f5c4e33-1c8b-41a1-ba73-52d3643fa008


# ╔═╡ Cell order:
# ╠═915f4df0-bc3a-11ec-2560-c9772e143679
# ╠═60d7decb-a7fa-494e-93ca-26d8b957cfcf
# ╠═7603c764-9e3f-4a66-9ee6-e7694e5070a3
# ╠═9175d027-b5eb-4676-bed3-e084a882e102
# ╠═aad0407a-4979-4ea4-94d3-38f1bb6d3de1
# ╠═943d7c37-ade5-4a4d-ae0b-e2bf390fcf6f
# ╠═85ddcf05-c9c5-4436-8511-3713dbbe7694
# ╠═59b234fe-0376-497c-bea3-997ab50422c5
# ╠═8493b4cb-12fa-4a1b-9125-ea59523a889b
# ╠═61a6b383-b9e2-4acf-8bb7-0fa926a0ec12
# ╠═55e42e79-b7f1-4793-aaa7-38ecfe0f3b66
# ╠═88313435-2633-40e3-94b1-661c15be2b30
# ╠═a984a49c-8fca-435c-a49b-09e72d5f6c48
# ╠═c7a089e6-9cb1-4a18-a65d-ff4044204eec
# ╠═ba40741a-3c3f-4e22-9ad1-4c150db5deb8
# ╠═8058af4a-4df7-45e3-82ff-53bf10908baf
# ╠═1e2aeb26-5f62-4f06-956a-1926befb923a
# ╠═b9c09401-186e-43c3-98db-f36538a88616
# ╠═bcd3233e-2e0a-44f0-8f18-9123a9a89012
# ╠═323a8caf-bb26-4c7b-b6b3-524351a80e5c
# ╠═6e7cca77-de89-42c5-a8ce-c016eb327fe3
# ╟─cfc6c757-8893-4d9c-add8-4fd6906ae156
# ╟─09444324-aa7b-4e63-af94-fa8950ece8b6
# ╟─287a44bb-e4d9-4a8b-b383-13d9c01e6e1e
# ╟─f0f007c1-6528-4a85-a007-b31881297f6a
# ╠═c7a6ed94-e97a-4a96-8f07-1370e71038b9
# ╠═e2c11399-7878-4602-a276-e190857b3fa6
# ╠═6dec59bf-1fa0-4ad6-afdd-26ec2b92230a
# ╠═560fda7a-4f4e-49e9-bf09-d82557286a83
# ╟─2c24347b-9c55-4e14-876a-ea8e2867eaa2
# ╠═f0767520-04f6-4e97-a889-cf5e45d70b4c
# ╠═a3eca5f4-21c0-4147-b7d0-19ba3fd7a776
# ╟─12a51952-34ab-47a3-a318-c4f58e3f3c12
# ╠═2cc7237d-3a8b-4b86-bb08-6cba2e9e9531
# ╠═5b66805c-42ef-4611-a31b-962ba3f549ed
# ╠═4653e385-beb8-49d6-a562-aae001db013e
# ╟─c6720b8a-e709-4c47-944f-dd4e7d750cbb
# ╠═454f4cb1-1535-4aba-a824-181489d324de
# ╠═2f2d6daf-c0b8-4501-86a5-65a844c9359e
# ╟─e900b657-4898-4811-aa6c-e178d693e6bf
# ╠═5aebef7b-6041-48b0-90cc-155975667b34
# ╟─bc5d6541-af3e-4101-a18e-fdcb8f3a7871
# ╠═1408b638-47ec-4cb7-9cd2-a2fea531eb8f
# ╟─2f34272e-68d1-4356-ac7e-13aa4e221b7d
# ╠═c04f6eaa-3c8a-4c8e-9311-63eb2e259315
# ╟─1d940676-91a4-4c58-a390-9fc276d706bc
# ╟─d711e132-192c-47c1-b8a0-7070aa5f55f4
# ╟─a6e855e1-9992-4969-b56e-872866801a9b
# ╟─be1f1295-a9d5-4dc8-9b6b-8269bbc15827
# ╟─03bf51d8-9735-458e-9d3a-816c9fcdde0c
# ╟─91ad223d-c245-4aca-89ab-cd6ee8c6652d
# ╟─d0b5715e-1964-4957-884b-5315f2a61777
# ╟─cb71ea78-9dc1-4bf9-871a-1a7f503f9371
# ╟─7b588922-125b-4866-868d-3ba9e63a307b
# ╠═26f3d3ee-91ab-42f1-bc94-c8a7c689e1d2
# ╟─e2325e9e-94bc-490f-877b-0922757fbd62
# ╟─56bbf4fd-c7ff-4532-bfd3-12f802a9bbd5
# ╟─923fc804-c188-417e-9465-2380b139469b
# ╠═45811f40-ed90-4155-938a-f6fd91f7030c
# ╟─da4e2714-7aec-4838-81a9-621c403130fb
# ╠═ec05e159-50de-439a-8063-a9add8db5599
# ╟─37a5cdeb-f98d-4eb0-9e5b-7e0935a26316
# ╠═1d3f03de-5b5b-4ea9-839c-a6ae67f72a22
# ╠═35a87a6e-4c70-400f-9e94-7b48e90e536f
# ╟─7e35a352-bc77-4eb1-9e5b-1a87fba90cc1
# ╟─efa579a9-f8ee-4158-9ed5-1c17295a47c9
# ╟─a605d131-3e89-4701-a660-f32ad6182811
# ╟─39e2a785-e432-40a8-b679-49e9440cbf99
# ╠═8956e904-e17f-4f9c-876f-8783dbd48bcd
# ╟─fc3b4271-4ea8-4e1f-a62f-69186d4fbac8
# ╟─50ba48db-8466-4299-967e-825ab99fba9e
# ╠═dcaf7bd1-c692-4067-a57c-d32f244a3054
# ╟─d65328b2-34f1-4321-b7d2-0dabca98e2e2
# ╠═4d0e2692-cd10-455b-9030-7b59046504a9
# ╠═7f5c4e33-1c8b-41a1-ba73-52d3643fa008
