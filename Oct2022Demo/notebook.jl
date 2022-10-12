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
include("Oct2022Demo.jl");

# ╔═╡ 287a44bb-e4d9-4a8b-b383-13d9c01e6e1e
md"""## Model Stratification"""

# ╔═╡ 4330976d-ceeb-4a75-97d0-a8375ede795b
md"""### Load Disease Model from JSON"""

# ╔═╡ 80e375d8-60fc-4857-b430-9973117c5d29
# mdl_disease_w_rates = read_json_acset(LabelledReactionNet{Any,Any},"../data/SIR.json");

# ╔═╡ 67264e96-d893-4f97-b433-ec24cd14d8ef
# mdl_disease = LabelledPetriNet(mdl_disease_w_rates);

# ╔═╡ b650ef6c-7138-4171-9481-288654d9a480
md"""**Form ACSet Type to Load MIRA JSON**"""

# ╔═╡ 9806cc7e-f14b-453f-af9a-0c0c87e6559c
begin
	@present TheoryOrigMIRANet <: TheoryLabelledReactionNet begin
    	MID::AttrType
    	MCTX::AttrType
    	Template::AttrType
    	mira_ids::Attr(S, MID)
    	mira_context::Attr(S, MCTX)
    	template_type::Attr(T, Template)
    	parameter_name::Attr(T, Name)
    	parameter_value::Attr(T, Rate)
	end
	@abstract_acset_type AbstractOrigMIRANet <: AbstractLabelledReactionNet
	@acset_type OrigMIRANet(TheoryOrigMIRANet) <: AbstractOrigMIRANet
end;

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

# ╔═╡ 02a79d90-34dc-4d6c-b2c8-388acfc5ab95
md"""**Load and Graph Original MIRA JSON**"""

# ╔═╡ 562f06ac-e2df-4040-942c-5d8494dbab4f
mdl_orig_mira = read_json_acset(OrigMIRANet{Any,Any,Any,Any,Any,Any},"BIOMD0000000971_petri_orig.json");

# ╔═╡ 173f8103-edc8-44fe-a2f0-79bfbdcd78d5
begin
	mdl_orig_mira[:sname] .= Symbol.(mdl_orig_mira[:sname])
	mdl_orig_mira[:tname] .= Symbol.(mdl_orig_mira[:tname])
	mdl_orig = LabelledPetriNet(mdl_orig_mira);
end;


# ╔═╡ 455cbff2-c001-4427-89ee-87393d978d06
AlgebraicPetri.Graph(mdl_orig)

# ╔═╡ f0f007c1-6528-4a85-a007-b31881297f6a
md"""**Load Disease Model from Curated MIRA JSON**"""

# ╔═╡ e2c11399-7878-4602-a276-e190857b3fa6
mdl_disease_mira = read_json_acset(MIRANet{Any,Any,Any,Any,Any,Any},"BIOMD0000000971_petri_curated.json");

# ╔═╡ a523d826-3d49-41eb-af52-849d24f0d919
md"""**Fill Rates and Concentrations**"""

# ╔═╡ 92ba176d-e991-41fd-88bf-191645ed9c88
begin
	input_rates = deepcopy(mdl_disease_mira[:rate])
	min_rate = minimum(input_rates[map(!isnothing,input_rates)])
	filled_rates = input_rates
	filled_rates[map(isnothing,input_rates)] = repeat([0.1],sum(map(isnothing,input_rates)))
end;

# ╔═╡ 72eb264f-061b-4360-a5bd-a12633f33a34
begin
	mdl_disease_mira[:rate] = filled_rates
	mdl_disease_mira[:concentration] = 0.0
	mdl_disease_mira[:sname] .= Symbol.(mdl_disease_mira[:sname])
	mdl_disease_mira[:tname] .= Symbol.(mdl_disease_mira[:tname])
end;

# ╔═╡ a2369b71-43a6-42bf-b322-82e56da8c151
md"""**Form Labelled Reaction Net**"""

# ╔═╡ 6dec59bf-1fa0-4ad6-afdd-26ec2b92230a
begin
	mdl_disease_rxnet = LabelledReactionNet{Any,Any}()
	copy_parts!(mdl_disease_rxnet,mdl_disease_mira)
	mdl_disease_rxnet
end

# ╔═╡ 7c16526d-7d9b-466a-b8de-6dc2e9c7a74f
mdl_disease_rxnet[:sname]

# ╔═╡ 55cb5f05-2613-48a7-b72f-9d4c94eccaee
md"""**Form Labelled Petri Net**"""

# ╔═╡ 560fda7a-4f4e-49e9-bf09-d82557286a83
mdl_disease = LabelledPetriNet(mdl_disease_mira);

# ╔═╡ 0e104aa0-07c8-4870-a976-7fc6cf8c25de
AlgebraicPetri.Graph(mdl_disease)

# ╔═╡ 2c24347b-9c55-4e14-876a-ea8e2867eaa2
md"""**Define Type System**"""

# ╔═╡ f0767520-04f6-4e97-a889-cf5e45d70b4c
begin
	types′ = LabelledPetriNet([:Pop],
    	:infect=>((:Pop, :Pop)=>(:Pop, :Pop)),
    	:disease=>(:Pop=>:Pop),
    	:strata=>(:Pop=>:Pop))
	types = map(types′, Name=name->nothing)
end;

# ╔═╡ 7c003649-03c0-4793-a673-d380e9d199ae
AlgebraicPetri.Graph(types′)

# ╔═╡ fbe5adf1-1d6b-40a6-b36e-49e3fc76057e
md"""**Define Typed Disease Model**"""

# ╔═╡ 4c3d24e2-254e-427e-9eb5-2344c1e8406d
#=begin 
	mdl_disease_typed = homomorphism(mdl_disease, types;
    	initial=(T=[1,2],I=[1,2,3],O=[1,2,3]),
    	type_components=(Name=x->nothing,))
	@assert is_natural(mdl_disease_typed)
end=#

# ╔═╡ a3eca5f4-21c0-4147-b7d0-19ba3fd7a776
begin 
	mdl_disease_typed = homomorphism(mdl_disease, types;
    	initial=(T=[1,1,2,2,2,2,2,2,2,2,2,1,2,1,1,1],),
    	type_components=(Name=x->nothing,))
	@assert is_natural(mdl_disease_typed)
end

# ╔═╡ 12a51952-34ab-47a3-a318-c4f58e3f3c12
md"""### Draw Stratification Model"""

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

# ╔═╡ 5b66805c-42ef-4611-a31b-962ba3f549ed
begin
	mdl_strat_typed = homomorphism(mdl_strat, types;
    	initial=(T=[3,3],), type_components=(Name=x->nothing,))
	@assert is_natural(mdl_strat_typed)
end;

# ╔═╡ b598e5b8-38e2-4d7d-afc9-2b0e7b1d2be8
md"""### Compute Stratified Model"""

# ╔═╡ ec24996f-12f8-4d82-b210-8b23529fa73b
md"""**Specify Stratification Structure**"""

# ╔═╡ be5b29e0-acb5-4bdd-951e-26fb1150a827
begin
	disease_ss = StrataSpec(mdl_disease_typed, [[:strata],[:strata],[:strata],[]])
	strat_ss = StrataSpec(mdl_strat_typed, [[:disease,:infect], [:disease,:infect]])
end;

# ╔═╡ 1a35abc0-66cf-41ae-b313-7dd248cef669
md"""**Stratify**"""

# ╔═╡ 2f2d6daf-c0b8-4501-86a5-65a844c9359e
mdl_stratified, obs_stratified = stratify(disease_ss, strat_ss, types′);

# ╔═╡ dad0382a-3656-4bc9-acf6-b914916b15ac
AlgebraicPetri.Graph(mdl_stratified)

# ╔═╡ 1c42f5e2-4e5c-4616-ab2d-c4d3ac877b9b
md"""### Save Stratified Model"""

# ╔═╡ cbcc0224-dae5-49c3-81eb-b1d37abbd526
# write_json_acset(mdl_stratified, "tmp_stratified_model.json")

# ╔═╡ 6b2b41bd-ea99-4917-8eec-d6963243f35a
md"""### Three Stratified Models"""

# ╔═╡ 1a5c2273-a1d1-44a4-a312-db8e9cd3c7de
md"""**SIRD Model**"""

# ╔═╡ 44328cee-fbf9-42db-afeb-7fc48815952e
SIRD = LabelledPetriNet([:S, :I, :R, :D],
    :inf => ((:S,:I) => (:I,:I)),
    :recover => (:I=>:R),
    :death => (:I=>:D));

# ╔═╡ 84fec776-b42c-4066-ad45-624b1eb93460
AlgebraicPetri.Graph(SIRD)

# ╔═╡ d482440e-2cce-4f3d-b4cb-208861a71187
begin 
	SIRD_typed = homomorphism(SIRD, types;
    	initial=(T=[1,2,2],I=[1,2,3,3],O=[1,2,3,3]),
    	type_components=(Name=x->nothing,))
	@assert is_natural(SIRD_typed)
end

# ╔═╡ f7d18d80-f2a6-474c-8115-f323e09637d8
md"""**SIRD-Quarantine Model**"""

# ╔═╡ c4f93cf2-9b30-44b2-9b7a-beee37686f37
begin
	Quarantine = LabelledPetriNet([:Q,:NQ],
    	:quarantine => ((:NQ)=>(:Q)),
    	:unquarantine => ((:Q)=>(:NQ)));
	Quarantine_typed = homomorphism(Quarantine, types;
    	initial=(T=[3,3],), type_components=(Name=x->nothing,))
	@assert is_natural(Quarantine_typed)
end

# ╔═╡ 3c922e81-61af-4789-b2f4-e96dd981810f
AlgebraicPetri.Graph(Quarantine)

# ╔═╡ cc4d55cc-3ad1-47ed-9755-2ef3c616bf57
begin
	SIRD_Q_ss = StrataSpec(SIRD_typed, [[:strata],[:strata],[:strata],[]])
	Quarantine_ss = StrataSpec(Quarantine_typed, [[:disease], [:disease,:infect]])
	mdl_SIRD_Q, obsSIRD_Q = stratify(SIRD_Q_ss, Quarantine_ss, types′);
end;

# ╔═╡ 0de945ad-093d-4975-8e48-a0748e2414b1
AlgebraicPetri.Graph(mdl_SIRD_Q)

# ╔═╡ 385a2551-4733-4fb0-9d5d-9763f6757578
md"""**SIRD-Mask Model (v1 - Kris's Stratification Form)**"""

# ╔═╡ 452ade8f-6140-4531-9919-51e899b15f33
begin
	Mask = LabelledPetriNet([:M,:NM],
    	:mask => ((:NM)=>(:M)),
    	:unmask => ((:M)=>(:NM)),
		:mu_infect => ((:M,:NM)=>(:M,:NM)),
		:uu_infect => ((:NM,:NM)=>(:NM,:NM)));
	Mask_typed = homomorphism(Mask, types;
    	initial=(T=[3,3,1,1],), type_components=(Name=x->nothing,))
	@assert is_natural(Mask_typed)
end

# ╔═╡ e79d741f-702d-49c9-9e1a-6442cfeb7614
AlgebraicPetri.Graph(Mask)

# ╔═╡ c33580a9-1d22-42a5-b527-018c2448b209
begin
	SIRD_M_ss = StrataSpec(SIRD_typed, [[:strata],[:strata],[:strata],[]])
	Mask_ss = StrataSpec(Mask_typed, [[:disease,:infect], [:disease,:infect]])
	mdl_SIRD_M, obsSIRD_M = stratify(SIRD_M_ss, Mask_ss, types′);
end;

# ╔═╡ a9043c9e-ec1b-41e8-a617-5a1a4d611415
AlgebraicPetri.Graph(mdl_SIRD_M)

# ╔═╡ 803517ca-ad26-4c3c-b003-78d8b9bdf0cc
md"""**SIRD-Mask Model (v2 - Sophie's Stratification Form)**"""

# ╔═╡ fb274f1d-fab3-4d50-ba87-c8f1cd880822
begin 
	s, = parts(types′, :S)
	t_interact, t_disease, t_strata = parts(types′, :T)
	i_interact1, i_interact2, i_disease, i_strata = parts(types′, :I)
	o_interact1, o_interact2, o_disease, o_strata = parts(types′, :O);
end;

# ╔═╡ 75d3e1cd-b951-4626-9625-9ef1c546ef87
begin 
	SIRD_aug = LabelledPetriNet([:S, :I, :R, :D],
	  :inf => ((:S, :I)=>(:I, :I)),
	  :rec => (:I=>:R),
	  :death => (:I=>:D),
	  :id => (:S => :S),
	  :id => (:I => :I),
	  :id => (:R => :R)
	#  :id => (:D => :D)
	)
	
	SIRD_aug_typed = ACSetTransformation(SIRD_aug, types,
	  S = [s, s, s, s],
	  T = [t_interact, t_disease, t_disease, t_strata, t_strata, t_strata],
	  I = [i_interact1, i_interact2, i_disease, i_disease, i_strata, i_strata, i_strata],
	  O = [o_interact1, o_interact2, o_disease, o_disease, o_strata, o_strata, o_strata],
	  Name = name -> nothing # specify the mapping for the loose ACSet transform
	)
	
	@assert is_natural(SIRD_aug_typed)
end;

# ╔═╡ c0458439-2af8-4e42-90c7-41f562d0750e
AlgebraicPetri.Graph(dom(SIRD_aug_typed)) # Graph_typed

# ╔═╡ 488e5c05-39ee-4635-9346-730bea5aa5f5
begin
	Mask_aug = LabelledPetriNet([:M,:NM],
    	:mu_infect => ((:M,:NM)=>(:M,:NM)),
		:uu_infect => ((:NM,:NM)=>(:NM,:NM)),
	  	:mask => ((:NM)=>(:M)),
    	:unmask => ((:M)=>(:NM)),
		:id => (:M => :M),
	  	:id => (:NM => :NM)
	)
	
	Mask_aug_typed = ACSetTransformation(Mask_aug, types,
	  S = [s, s],
	  T = [t_interact, t_interact, t_strata, t_strata, t_disease, t_disease],
	  I = [i_interact1, i_interact2, i_interact1, i_interact2, i_strata, i_strata, i_disease, i_disease],
	  O = [o_interact1, o_interact2, o_interact1, o_interact2, o_strata, o_strata, o_disease, o_disease],
	  Name = name -> nothing # specify the mapping for the loose ACSet transform
	)

	@assert is_natural(Mask_aug_typed)
end;

# ╔═╡ ece17ed4-dc74-4bfc-b440-6f34cd560cab
AlgebraicPetri.Graph(dom(Mask_aug_typed)) # Graph_typed

# ╔═╡ 06636e60-e7dc-43ee-8541-190881085d16
begin
	typed_stratify(typed_model1, typed_model2) =
		Theories.compose(proj1(CategoricalAlgebra.pullback(typed_model1, typed_model2)), typed_model1)
	mdl2_SIRD_M = typed_stratify(SIRD_aug_typed, Mask_aug_typed)
end;

# ╔═╡ 03ad5137-c6f7-424d-99fb-65ac7f8a2ed2
AlgebraicPetri.Graph(dom(mdl2_SIRD_M)) # Graph_typed

# ╔═╡ d9fd32c2-ed42-4af6-90a7-d144296e222d
md"""**SIRD-Two-CityModel**"""

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

# ╔═╡ 38f0c179-16b2-4c4d-9975-8fb4c3f9f97b
begin
	SIRD_TC_ss = StrataSpec(SIRD_typed, [[:strata],[:strata],[:strata],[]])
	TwoCity_ss = StrataSpec(TwoCity_typed, [[:disease,:infect], [:disease,:infect]])
	mdl_SIRD_TC, obsSIRD_TC = stratify(SIRD_TC_ss, TwoCity_ss, types′);
end;

# ╔═╡ 015fa0fd-5e9f-469a-bbb8-c388b946b3ab
AlgebraicPetri.Graph(mdl_SIRD_TC)

# ╔═╡ cfc6c757-8893-4d9c-add8-4fd6906ae156
md"""## Model Calibration"""

# ╔═╡ cabe7fdb-73bd-493b-80d2-0fc8f1f38be7
md"""### Generate Data"""

# ╔═╡ 157be78c-c73f-484b-b944-a733f811457a
#=begin
	true_mdl = mdl_SIRD_Q
	true_obs = obsSIRD_Q
	true_p = [5.0e-4, 5.0e-4, 5.0e-4, 5.0e-4, 5.0e-4, 5.0e-4,
            1.0e-2, 1.0e-5, 1.0e-3, 1.0e-5, 5.0e-5]
	u0 = [0.0,0.0,0.0,0.0,999.0,1.0,0.0,0.0] 
	tspan = (0.0,250.0)
	sample_data, sample_times, prob_real, true_sol, noiseless_data, data_labels = 		generateData(true_mdl, true_p, u0, tspan, 50, true_obs)
end;=#

# ╔═╡ 13f557c7-9002-4a2d-8be5-fccef5aa9dc6
md"""### Calibrate"""

# ╔═╡ 3b8954ae-2959-4ced-bdd2-1381121eb6c9
#=begin
	states_to_count = [:S,:I,:R,:D]
	p_init = repeat([1.0e-6], nt(true_mdl))
	p_est, sol_est, loss = calibrate(true_mdl, true_obs, states_to_count, u0, p_init, 	sample_data, sample_times, data_labels)
end;=#

# ╔═╡ 2276bfa9-c6bd-4c54-a6e1-7a88fd7a7762
md"""### Results"""

# ╔═╡ 79633d27-71b4-4261-81dc-db609133e2cd
md"""**Estimated Parameters**"""

# ╔═╡ ef69c329-9e2d-444a-8fa5-447f15b419b1
# hcat(true_mdl[:tname],p_est)

# ╔═╡ 60408568-5762-43aa-9fef-b10a8c903422
md"""**Estimated Observations**"""

# ╔═╡ 8f1bfc19-0c0d-4ac9-a16a-fb56555f6707
# plot_obs_w_ests(sample_times, sample_data, sol_est, true_obs)

# ╔═╡ 1d940676-91a4-4c58-a390-9fc276d706bc
md"""## Elaborate Model"""

# ╔═╡ d512194b-6811-4433-9383-ccc9d78463f9
#=Elaborate = LabelledPetriNet([:S, :E, :I, :A, :SQ, :H, :R, :EQ, :D],
	:expos_a => ((:S,:A) => (:E,:A)),
	:spook_sa => ((:S,:A)=>(:SQ,:A)),
	:unspook_s => (:SQ=>:S),
	:prog_ei => (:E=>:I),
	:prog_ea => (:E=>:A),
	:hosp_i => (:I=>:H),
	:recov_i => (:I=>:R),
	:recov_a => (:A=>:R),
	:recov_h => (:H=>:R),
	:death_i => (:I=>:D),
	:death_h => (:H=>:D),	
	:espook_a => ((:S,:A)=>(:EQ,:A)),
	:hosp_eq => (:EQ=>:H),
	#:spook_ea => ((:E,:A)=>(:EQ,:A)),
	#:spook_ei => ((:E,:I)=>(:EQ,:I)),	
	:expos_i => ((:S,:I) => (:E,:I)),
	:spook_si => ((:S,:I)=>(:SQ,:I)),
	:espook_i => ((:S,:I)=>(:EQ,:I))
);=#

# ╔═╡ 8f985d35-d601-43d8-8e66-aa035f8ee63d
Elaborate = mdl_disease;

# ╔═╡ 1408b638-47ec-4cb7-9cd2-a2fea531eb8f
AlgebraicPetri.Graph(Elaborate)

# ╔═╡ d711e132-192c-47c1-b8a0-7070aa5f55f4
md"""**Define Observation Function**"""

# ╔═╡ 394b9c92-3681-4064-b5a7-a62d83814dc7
begin 
function stateidx(model, name)
	filter(parts(model, :S)) do i
		model[i, :sname] == name
	end
end

function stateidx_stratified(model, name, dim=1)
	filter(parts(model, :S)) do i
		model[i, :sname][dim] == name
	end
end

function sumvarsbyname(model, name, sol, sample_times)
	idxs = statidx(model, :I)
    sample_vals = sum(sol(sample_times)[idxs,:], dims=1)
end
end;

# ╔═╡ 234195c0-ebd8-4214-af62-5c486006a814
function obs_Elaborate(model::AbstractLabelledPetriNet, sol, sample_times)
    inf_sample_vals = sumvarsbyname(model, :I, sol, sample_times)
    hosp_sample_vals = sumvarsbyname(model, :H, sol, sample_times)
	dead_sample_vals = sumvarsbyname(model, :D, sol, sample_times)
	
    labels = reshape(["I", "H", "D"],1,3)

    return hcat(inf_sample_vals, hosp_sample_vals, dead_sample_vals), labels
end;

# ╔═╡ 01b6524d-6b68-4dd2-8cc5-ffaf9678835a
#=begin 
	Elaborate_typed = homomorphism(Elaborate, types;
    	initial=(T=[1,1,2,2,2,2,2,2,2,2,2,1,2,1,1,1],),
    	type_components=(Name=x->nothing,))
	@assert is_natural(Elaborate_typed)
end=#

# ╔═╡ c82f017c-0da4-442a-a25e-36fc5c7debfa
Elaborate_typed = mdl_disease_typed;

# ╔═╡ 973aadbd-9f19-4e17-9b7a-bb977b405d4c
#=begin
	Elaborate_TC_ss = StrataSpec(Elaborate_typed, 
		[[:strata],[:strata],[:strata],[:strata],[],[],[:strata],[],[]])
	TwoCity_Elab_ss = StrataSpec(TwoCity_typed, [[:disease,:infect], 		[:disease,:infect]])
	mdl_Elab_TC, obs_Elab_TC = stratify(Elaborate_TC_ss, TwoCity_Elab_ss, types′);
end;=#

# ╔═╡ 5aebef7b-6041-48b0-90cc-155975667b34
begin
	mdl_Elab_TC = mdl_stratified 
	obs_Elab_TC = obs_stratified
end;

# ╔═╡ c04f6eaa-3c8a-4c8e-9311-63eb2e259315
AlgebraicPetri.Graph(mdl_Elab_TC)

# ╔═╡ 26f3d3ee-91ab-42f1-bc94-c8a7c689e1d2
function obs_Elab_TC_trunc(obs_func, sol, sample_times)
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

# ╔═╡ d9c4102f-f8de-46e3-ab40-4a8e893a7a0c
#=function obs_Elab_TC_trunc_old(diseasemodel::AbstractLabelledPetriNet, sol, sample_times)
	samples, labels = obs_Elab_TC(sol, sample_times)
	indx(labels, l) = findall(isequal(l), labels)
	sumsamples(l) = sum(samples[:, indx(labels, l)], dims=2)
	i = sumsamples(:I)
	h = sumsamples(:H)
	d = sumsamples(:D)
	labels = [:I :H :D]
	# labels = reshape(["I", "H", "D"],1,3)

    return hcat(i, h, d), labels
end;=#

# ╔═╡ 45811f40-ed90-4155-938a-f6fd91f7030c
begin
	true_mdl = mdl_Elab_TC
	# true_obs = (sol, times) -> obs_Elab_TC_trunc(obs_Elab_TC, sol, times)
	true_obs = obs_Elab_TC
	u0 = repeat([999.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0],2)
end;

# ╔═╡ ec05e159-50de-439a-8063-a9add8db5599
begin 
	function transidx(model, name)
		filter(parts(model, :T)) do i
			model[i, :tname] == name
		end
	end
	true_p = zeros(nt(true_mdl))
	for ii in 1:nt(true_mdl)
		if true_mdl[ii,:tname][1]==:strata
			true_p[ii] = .2 # true_mdl[ii,:tname][1]
		else 
			true_p[ii] = filled_rates[transidx(Elaborate,true_mdl[ii,:tname][1])[1]]
		end
	end
end;

# ╔═╡ 1d3f03de-5b5b-4ea9-839c-a6ae67f72a22
begin
	#true_p = repeat(filled_rates,2)
	tspan = (0.0,250.0)
	sample_data, sample_times, prob_real, true_sol, noiseless_data, data_labels = 		generateData(true_mdl, true_p, u0, tspan, 50, true_obs)
end;

# ╔═╡ 66d94ce5-9ca3-4597-95ca-acc0722495eb
hcat(Elaborate[:tname],filled_rates)

# ╔═╡ 5c896c1e-e229-45c2-b2db-f81ae5d80fa5
sample_data |> size

# ╔═╡ b868463a-dc90-4fa1-b3ae-bf8ff11824be
maximum(sample_data[:,3])

# ╔═╡ 35a87a6e-4c70-400f-9e94-7b48e90e536f
plot(sample_data, labels=String.(data_labels))

# ╔═╡ 4cdc078f-ee7a-4c05-aaeb-041ac4b85f74
true_mdl[:tname]

# ╔═╡ 8956e904-e17f-4f9c-876f-8783dbd48bcd
begin
	states_to_count = [:I,:H,:D]
	# p_init = repeat([1.0e-6], nt(true_mdl))
	p_init = true_p
	p_est, sol_est, loss = calibrate(true_mdl, true_obs, states_to_count, u0, p_init, 	sample_data, sample_times, data_labels)
end;

# ╔═╡ dcaf7bd1-c692-4067-a57c-d32f244a3054
hcat(true_mdl[:tname],p_est)

# ╔═╡ 4d0e2692-cd10-455b-9030-7b59046504a9
plot_obs_w_ests(sample_times, sample_data, sol_est, true_obs)

# ╔═╡ f7664021-998c-403b-bf73-c5bec585854f
apex[:rate] .= *(apex[:rate][1], apex[:rate][2])

# ╔═╡ Cell order:
# ╠═915f4df0-bc3a-11ec-2560-c9772e143679
# ╠═60d7decb-a7fa-494e-93ca-26d8b957cfcf
# ╠═aad0407a-4979-4ea4-94d3-38f1bb6d3de1
# ╟─287a44bb-e4d9-4a8b-b383-13d9c01e6e1e
# ╟─4330976d-ceeb-4a75-97d0-a8375ede795b
# ╟─80e375d8-60fc-4857-b430-9973117c5d29
# ╟─67264e96-d893-4f97-b433-ec24cd14d8ef
# ╟─b650ef6c-7138-4171-9481-288654d9a480
# ╠═9806cc7e-f14b-453f-af9a-0c0c87e6559c
# ╠═c7a6ed94-e97a-4a96-8f07-1370e71038b9
# ╟─02a79d90-34dc-4d6c-b2c8-388acfc5ab95
# ╠═562f06ac-e2df-4040-942c-5d8494dbab4f
# ╠═173f8103-edc8-44fe-a2f0-79bfbdcd78d5
# ╠═455cbff2-c001-4427-89ee-87393d978d06
# ╟─f0f007c1-6528-4a85-a007-b31881297f6a
# ╠═e2c11399-7878-4602-a276-e190857b3fa6
# ╟─a523d826-3d49-41eb-af52-849d24f0d919
# ╠═92ba176d-e991-41fd-88bf-191645ed9c88
# ╠═72eb264f-061b-4360-a5bd-a12633f33a34
# ╟─a2369b71-43a6-42bf-b322-82e56da8c151
# ╠═6dec59bf-1fa0-4ad6-afdd-26ec2b92230a
# ╠═7c16526d-7d9b-466a-b8de-6dc2e9c7a74f
# ╟─55cb5f05-2613-48a7-b72f-9d4c94eccaee
# ╠═560fda7a-4f4e-49e9-bf09-d82557286a83
# ╠═0e104aa0-07c8-4870-a976-7fc6cf8c25de
# ╟─2c24347b-9c55-4e14-876a-ea8e2867eaa2
# ╠═f0767520-04f6-4e97-a889-cf5e45d70b4c
# ╠═7c003649-03c0-4793-a673-d380e9d199ae
# ╟─fbe5adf1-1d6b-40a6-b36e-49e3fc76057e
# ╟─4c3d24e2-254e-427e-9eb5-2344c1e8406d
# ╠═a3eca5f4-21c0-4147-b7d0-19ba3fd7a776
# ╟─12a51952-34ab-47a3-a318-c4f58e3f3c12
# ╠═2cc7237d-3a8b-4b86-bb08-6cba2e9e9531
# ╠═597f6c5d-c404-4288-baf1-c892aa59deb2
# ╟─52cc571c-7d29-4167-834c-6672bbef6edf
# ╠═5b66805c-42ef-4611-a31b-962ba3f549ed
# ╟─b598e5b8-38e2-4d7d-afc9-2b0e7b1d2be8
# ╟─ec24996f-12f8-4d82-b210-8b23529fa73b
# ╠═be5b29e0-acb5-4bdd-951e-26fb1150a827
# ╟─1a35abc0-66cf-41ae-b313-7dd248cef669
# ╠═2f2d6daf-c0b8-4501-86a5-65a844c9359e
# ╠═dad0382a-3656-4bc9-acf6-b914916b15ac
# ╟─1c42f5e2-4e5c-4616-ab2d-c4d3ac877b9b
# ╠═cbcc0224-dae5-49c3-81eb-b1d37abbd526
# ╟─6b2b41bd-ea99-4917-8eec-d6963243f35a
# ╟─1a5c2273-a1d1-44a4-a312-db8e9cd3c7de
# ╠═44328cee-fbf9-42db-afeb-7fc48815952e
# ╠═84fec776-b42c-4066-ad45-624b1eb93460
# ╠═d482440e-2cce-4f3d-b4cb-208861a71187
# ╟─f7d18d80-f2a6-474c-8115-f323e09637d8
# ╠═c4f93cf2-9b30-44b2-9b7a-beee37686f37
# ╠═3c922e81-61af-4789-b2f4-e96dd981810f
# ╠═cc4d55cc-3ad1-47ed-9755-2ef3c616bf57
# ╠═0de945ad-093d-4975-8e48-a0748e2414b1
# ╟─385a2551-4733-4fb0-9d5d-9763f6757578
# ╠═452ade8f-6140-4531-9919-51e899b15f33
# ╠═e79d741f-702d-49c9-9e1a-6442cfeb7614
# ╠═c33580a9-1d22-42a5-b527-018c2448b209
# ╠═a9043c9e-ec1b-41e8-a617-5a1a4d611415
# ╟─803517ca-ad26-4c3c-b003-78d8b9bdf0cc
# ╠═fb274f1d-fab3-4d50-ba87-c8f1cd880822
# ╠═75d3e1cd-b951-4626-9625-9ef1c546ef87
# ╠═c0458439-2af8-4e42-90c7-41f562d0750e
# ╠═488e5c05-39ee-4635-9346-730bea5aa5f5
# ╠═ece17ed4-dc74-4bfc-b440-6f34cd560cab
# ╠═06636e60-e7dc-43ee-8541-190881085d16
# ╠═03ad5137-c6f7-424d-99fb-65ac7f8a2ed2
# ╟─d9fd32c2-ed42-4af6-90a7-d144296e222d
# ╠═4c7d9df5-8f92-4fe1-a2e3-82f607d28ccb
# ╠═25ba5ac6-9c41-4dfc-93b2-75664ea59565
# ╠═38f0c179-16b2-4c4d-9975-8fb4c3f9f97b
# ╠═015fa0fd-5e9f-469a-bbb8-c388b946b3ab
# ╟─cfc6c757-8893-4d9c-add8-4fd6906ae156
# ╟─cabe7fdb-73bd-493b-80d2-0fc8f1f38be7
# ╠═157be78c-c73f-484b-b944-a733f811457a
# ╟─13f557c7-9002-4a2d-8be5-fccef5aa9dc6
# ╠═3b8954ae-2959-4ced-bdd2-1381121eb6c9
# ╟─2276bfa9-c6bd-4c54-a6e1-7a88fd7a7762
# ╟─79633d27-71b4-4261-81dc-db609133e2cd
# ╠═ef69c329-9e2d-444a-8fa5-447f15b419b1
# ╟─60408568-5762-43aa-9fef-b10a8c903422
# ╠═8f1bfc19-0c0d-4ac9-a16a-fb56555f6707
# ╟─1d940676-91a4-4c58-a390-9fc276d706bc
# ╠═d512194b-6811-4433-9383-ccc9d78463f9
# ╠═8f985d35-d601-43d8-8e66-aa035f8ee63d
# ╠═1408b638-47ec-4cb7-9cd2-a2fea531eb8f
# ╟─d711e132-192c-47c1-b8a0-7070aa5f55f4
# ╠═394b9c92-3681-4064-b5a7-a62d83814dc7
# ╠═234195c0-ebd8-4214-af62-5c486006a814
# ╠═01b6524d-6b68-4dd2-8cc5-ffaf9678835a
# ╠═c82f017c-0da4-442a-a25e-36fc5c7debfa
# ╠═973aadbd-9f19-4e17-9b7a-bb977b405d4c
# ╠═5aebef7b-6041-48b0-90cc-155975667b34
# ╠═c04f6eaa-3c8a-4c8e-9311-63eb2e259315
# ╠═26f3d3ee-91ab-42f1-bc94-c8a7c689e1d2
# ╠═d9c4102f-f8de-46e3-ab40-4a8e893a7a0c
# ╠═45811f40-ed90-4155-938a-f6fd91f7030c
# ╠═ec05e159-50de-439a-8063-a9add8db5599
# ╠═1d3f03de-5b5b-4ea9-839c-a6ae67f72a22
# ╠═66d94ce5-9ca3-4597-95ca-acc0722495eb
# ╠═5c896c1e-e229-45c2-b2db-f81ae5d80fa5
# ╠═b868463a-dc90-4fa1-b3ae-bf8ff11824be
# ╠═35a87a6e-4c70-400f-9e94-7b48e90e536f
# ╠═4cdc078f-ee7a-4c05-aaeb-041ac4b85f74
# ╠═8956e904-e17f-4f9c-876f-8783dbd48bcd
# ╠═dcaf7bd1-c692-4067-a57c-d32f244a3054
# ╠═4d0e2692-cd10-455b-9030-7b59046504a9
# ╠═f7664021-998c-403b-bf73-c5bec585854f
