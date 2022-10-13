### A Pluto.jl notebook ###
# v0.19.13

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
include("Oct2022Demo.jl");

# ╔═╡ 287a44bb-e4d9-4a8b-b383-13d9c01e6e1e
md"""## Model Stratification"""

# ╔═╡ 4330976d-ceeb-4a75-97d0-a8375ede795b
md"""### Load Disease Model from JSON"""

# ╔═╡ 2cf45122-0246-4a07-9efb-e12c76fa9d98
md"""**Load Disease Model from Original MIRA JSON**"""

# ╔═╡ 9efe6732-77d5-4659-8ccd-4b01583454a3
md"""To read in the MIRA model from the JSON file we need to define a new ACSet Type."""

# ╔═╡ f0f007c1-6528-4a85-a007-b31881297f6a
md"""**Load Disease Model from Curated MIRA JSON**"""

# ╔═╡ 0f12d42f-1a6a-4ca0-bd60-3cb9f38fa518
md"""The original MIRA format is difficult for visualization, so we curated the JSON file. 
This curated format also required its own ACSet Type."""

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

# ╔═╡ a523d826-3d49-41eb-af52-849d24f0d919
md"""After reading in, rates and concentrations missing from the model still need to be filled. 
We provide default values."""

# ╔═╡ e2c11399-7878-4602-a276-e190857b3fa6
mdl_disease_mira = load_mira_curated("BIOMD0000000971_petri_curated.json", 0.1)

# ╔═╡ a2369b71-43a6-42bf-b322-82e56da8c151
md"""To perform the stratifications and other computations, we extract the model data into the LabelledReactionNet and LabelledPetriNet types, removing the extra MIRA fields."""

# ╔═╡ 6dec59bf-1fa0-4ad6-afdd-26ec2b92230a
begin
	filled_rates = mdl_disease_mira[:rate]
	mdl_disease_rxnet = LabelledReactionNet{Any,Any}()
	copy_parts!(mdl_disease_rxnet,mdl_disease_mira)
end;

# ╔═╡ 560fda7a-4f4e-49e9-bf09-d82557286a83
mdl_disease = LabelledPetriNet(mdl_disease_rxnet);

# ╔═╡ 2c24347b-9c55-4e14-876a-ea8e2867eaa2
md"""### Define Type System"""

# ╔═╡ 23f214d0-ab7c-4325-8ea5-25591f5d571c
md"""Performing stratification requires specifying the types of the components of the models to be stratified. This requires specification of a type system."""

# ╔═╡ b65b454d-8347-40e0-a4fd-44a353f3fe59
md"""We define a type system with one state and three transitions."""

# ╔═╡ f0767520-04f6-4e97-a889-cf5e45d70b4c
begin
	types′ = LabelledPetriNet([:Pop],
    	:infect=>((:Pop, :Pop)=>(:Pop, :Pop)),
    	:disease=>(:Pop=>:Pop),
    	:strata=>(:Pop=>:Pop))
	types = map(types′, Name=name->nothing)
end;

# ╔═╡ fbe5adf1-1d6b-40a6-b36e-49e3fc76057e
md"""**Define Typed Disease Model**"""

# ╔═╡ 8ddeb97d-eb39-4a9d-9eaf-fa955c9019e5
md"""We can now specify the typing of the MIRA disease model by specifying a homomorphism into the type system model."""

# ╔═╡ 4f0fc639-ae18-4c4f-9bd8-aa90037cb759
md"""In this case, the six interaction transitions are of infect type and all other transitions are of disease type."""

# ╔═╡ a3eca5f4-21c0-4147-b7d0-19ba3fd7a776
begin 
	mdl_disease_typed = homomorphism(mdl_disease, types;
    	initial=(T=[1,1,2,2,2,2,2,2,2,2,2,1,2,1,1,1],),
    	type_components=(Name=x->nothing,))
	@assert is_natural(mdl_disease_typed)
end

# ╔═╡ 12a51952-34ab-47a3-a318-c4f58e3f3c12
md"""### Draw Stratification Model"""

# ╔═╡ c89f40f2-3aa9-45c0-9cc8-0abf0a67564f
md"""Here we use `Semagrams.jl` to form the stratification model.""" 

# ╔═╡ 67fd498c-0460-4673-a2f7-ff45f824be3a
md"""As an example, we choose a Two-City movement model with states `:City1` and `:City2` and transitions `:travel12` and `:travel21`."""

# ╔═╡ 2cc7237d-3a8b-4b86-bb08-6cba2e9e9531
begin
	mdl_strat = LabelledPetriNet([:City1,:City2],
	    :travel12 => ((:City1)=>(:City2)),
    	:travel21 => ((:City2)=>(:City1)))
end;

# ╔═╡ 597f6c5d-c404-4288-baf1-c892aa59deb2
AlgebraicPetri.Graph(mdl_strat)

# ╔═╡ 52cc571c-7d29-4167-834c-6672bbef6edf
md"""**Define Typed Stratification Model**"""

# ╔═╡ f311a863-7580-49ca-8200-d6f053e5118e
md"""As with the disease model, we specify the typing of the stratification model. Here the two transitions are of type strata."""

# ╔═╡ 5b66805c-42ef-4611-a31b-962ba3f549ed
begin
	mdl_strat_typed = homomorphism(mdl_strat, types;
    	initial=(T=[3,3],), type_components=(Name=x->nothing,))
	@assert is_natural(mdl_strat_typed)
end;

# ╔═╡ b598e5b8-38e2-4d7d-afc9-2b0e7b1d2be8
md"""### Compute Stratified Model"""

# ╔═╡ 1185894b-063f-402f-96e2-8852430ac4d5
md"""To compute the stratified model, we must fully specify the available transitions and their types."""

# ╔═╡ 6b4584c6-dafa-466f-9879-299556a51277
md"""**Stratification Structs**"""

# ╔═╡ 9baf00c0-1b37-43f8-8a96-f75f45b33a77
md"""Here, we use a high-level, simplified form to make this specification. In this form, the desired transitions are indicated via stratification structs."""

# ╔═╡ 89199dc4-128b-49be-988c-bc69f5b44d8c
md"""Specifically, for each state of each model, we specify the transition types of the other model in which each state participates."""

# ╔═╡ 74e034c5-8f41-47e0-aea1-bdbe401fecd6
md"""In this case, the SQ, H, EQ, and D states cannot travel between cities, but all other states can. These transitions are of type strata. Both cities participate in disease and infect transitions."""

# ╔═╡ 454f4cb1-1535-4aba-a824-181489d324de
begin
	disease_ss = StrataSpec(mdl_disease_typed, 
		[[:strata],[:strata],[:strata],[:strata],[],[],[:strata],[],[]])
	strat_ss = StrataSpec(mdl_strat_typed, [[:disease,:infect], 			[:disease,:infect]])
end;

# ╔═╡ 1a35abc0-66cf-41ae-b313-7dd248cef669
md"""**Stratified Model**"""

# ╔═╡ 107561be-dcff-4613-bb78-54dee05ed625
md"""From the stratification structs and the type system, we can form the stratified model."""

# ╔═╡ 2f2d6daf-c0b8-4501-86a5-65a844c9359e
mdl_stratified, obs_stratified = stratify(disease_ss, strat_ss, types′);

# ╔═╡ 7836fa75-5ee2-410b-aa12-eded84fac180
md"""The resulting stratified model is"""

# ╔═╡ dad0382a-3656-4bc9-acf6-b914916b15ac
AlgebraicPetri.Graph(mdl_stratified)

# ╔═╡ 1663ff26-27b4-40e7-95fe-20b8ad5333f5
md"""Note that from the stratification we also obtain a particular observation function given by the projection of the pullback onto the disease model, i.e., the observations are the disease states summed across the corresponding stratified states."""


# ╔═╡ d9fd32c2-ed42-4af6-90a7-d144296e222d
md"""**SIRD-Two-City Model**"""

# ╔═╡ 44e9d001-e0df-433a-a05a-00c57e4c97c5
md"""Here we stratify SIRD with the model of movement between Two Cities."""

# ╔═╡ 0cca8978-d7b4-496b-8e48-3e5bf14a4028
md"""As with the MIRA disease model, this stratification can be specified using the previous high-level form, with the Two-City model given by"""

# ╔═╡ 4c7d9df5-8f92-4fe1-a2e3-82f607d28ccb
begin
	TwoCity = LabelledPetriNet([:City1,:City2],
	    :travel12 => ((:City1)=>(:City2)),
    	:travel21 => ((:City2)=>(:City1)))
	TwoCity_typed = homomorphism(TwoCity, types;
	    initial=(T=[3,3],), type_components=(Name=x->nothing,))
	@assert is_natural(TwoCity_typed)
end

# ╔═╡ 25ba5ac6-9c41-4dfc-93b2-75664ea59565
AlgebraicPetri.Graph(TwoCity)

# ╔═╡ 91ad223d-c245-4aca-89ab-cd6ee8c6652d
md"""Resulting in the following stratified model"""

# ╔═╡ cfc6c757-8893-4d9c-add8-4fd6906ae156
md"""## Model Calibration of MIRA-Two-City"""

# ╔═╡ 09444324-aa7b-4e63-af94-fa8950ece8b6
md"""Here we simulate data from the stratified model (formed from the MIRA disease model and Semagrams two-city stratification model) and fit the model to the data sample."""

# ╔═╡ e900b657-4898-4811-aa6c-e178d693e6bf
md"""For clarity we rename the disease and stratified models."""

# ╔═╡ 5aebef7b-6041-48b0-90cc-155975667b34
begin
	Mira = mdl_disease
	Mira_typed = mdl_disease_typed
	mdl_Mira_TC = mdl_stratified 
	obs_Mira_TC = obs_stratified
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
	u0 = repeat([999.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0],2)
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
	tspan = (0.0,250.0)
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
	# p_init = repeat([1.0e-6], nt(true_mdl))
	p_init = true_p
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

# ╔═╡ Cell order:
# ╠═915f4df0-bc3a-11ec-2560-c9772e143679
# ╠═60d7decb-a7fa-494e-93ca-26d8b957cfcf
# ╠═aad0407a-4979-4ea4-94d3-38f1bb6d3de1
# ╟─287a44bb-e4d9-4a8b-b383-13d9c01e6e1e
# ╟─4330976d-ceeb-4a75-97d0-a8375ede795b
# ╟─2cf45122-0246-4a07-9efb-e12c76fa9d98
# ╟─9efe6732-77d5-4659-8ccd-4b01583454a3
# ╟─f0f007c1-6528-4a85-a007-b31881297f6a
# ╟─0f12d42f-1a6a-4ca0-bd60-3cb9f38fa518
# ╠═c7a6ed94-e97a-4a96-8f07-1370e71038b9
# ╟─a523d826-3d49-41eb-af52-849d24f0d919
# ╠═e2c11399-7878-4602-a276-e190857b3fa6
# ╟─a2369b71-43a6-42bf-b322-82e56da8c151
# ╠═6dec59bf-1fa0-4ad6-afdd-26ec2b92230a
# ╠═560fda7a-4f4e-49e9-bf09-d82557286a83
# ╟─2c24347b-9c55-4e14-876a-ea8e2867eaa2
# ╟─23f214d0-ab7c-4325-8ea5-25591f5d571c
# ╟─b65b454d-8347-40e0-a4fd-44a353f3fe59
# ╠═f0767520-04f6-4e97-a889-cf5e45d70b4c
# ╟─fbe5adf1-1d6b-40a6-b36e-49e3fc76057e
# ╟─8ddeb97d-eb39-4a9d-9eaf-fa955c9019e5
# ╟─4f0fc639-ae18-4c4f-9bd8-aa90037cb759
# ╠═a3eca5f4-21c0-4147-b7d0-19ba3fd7a776
# ╟─12a51952-34ab-47a3-a318-c4f58e3f3c12
# ╟─c89f40f2-3aa9-45c0-9cc8-0abf0a67564f
# ╟─67fd498c-0460-4673-a2f7-ff45f824be3a
# ╠═2cc7237d-3a8b-4b86-bb08-6cba2e9e9531
# ╠═597f6c5d-c404-4288-baf1-c892aa59deb2
# ╟─52cc571c-7d29-4167-834c-6672bbef6edf
# ╟─f311a863-7580-49ca-8200-d6f053e5118e
# ╠═5b66805c-42ef-4611-a31b-962ba3f549ed
# ╟─b598e5b8-38e2-4d7d-afc9-2b0e7b1d2be8
# ╟─1185894b-063f-402f-96e2-8852430ac4d5
# ╟─6b4584c6-dafa-466f-9879-299556a51277
# ╟─9baf00c0-1b37-43f8-8a96-f75f45b33a77
# ╟─89199dc4-128b-49be-988c-bc69f5b44d8c
# ╟─74e034c5-8f41-47e0-aea1-bdbe401fecd6
# ╠═454f4cb1-1535-4aba-a824-181489d324de
# ╟─1a35abc0-66cf-41ae-b313-7dd248cef669
# ╟─107561be-dcff-4613-bb78-54dee05ed625
# ╠═2f2d6daf-c0b8-4501-86a5-65a844c9359e
# ╟─7836fa75-5ee2-410b-aa12-eded84fac180
# ╠═dad0382a-3656-4bc9-acf6-b914916b15ac
# ╟─1663ff26-27b4-40e7-95fe-20b8ad5333f5
# ╟─d9fd32c2-ed42-4af6-90a7-d144296e222d
# ╟─44e9d001-e0df-433a-a05a-00c57e4c97c5
# ╟─0cca8978-d7b4-496b-8e48-3e5bf14a4028
# ╠═4c7d9df5-8f92-4fe1-a2e3-82f607d28ccb
# ╠═25ba5ac6-9c41-4dfc-93b2-75664ea59565
# ╟─91ad223d-c245-4aca-89ab-cd6ee8c6652d
# ╟─cfc6c757-8893-4d9c-add8-4fd6906ae156
# ╟─09444324-aa7b-4e63-af94-fa8950ece8b6
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
