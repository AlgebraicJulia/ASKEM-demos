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
mdl_disease = read_json_acset(LabelledPetriNet,"../data/SIR.json");

# ╔═╡ e2c11399-7878-4602-a276-e190857b3fa6
# test_mdl = read_json_acset(LabelledReactionNet,"../data/SIR.json")

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

# ╔═╡ fbe5adf1-1d6b-40a6-b36e-49e3fc76057e
md"""**Define Typed Disease Model**"""

# ╔═╡ 4c3d24e2-254e-427e-9eb5-2344c1e8406d
begin 
	mdl_disease_typed = homomorphism(mdl_disease, types;
    	initial=(T=[1,2],I=[1,2,3],O=[1,2,3]),
    	type_components=(Name=x->nothing,))
	@assert is_natural(mdl_disease_typed)
end

# ╔═╡ 12a51952-34ab-47a3-a318-c4f58e3f3c12
md"""### Draw Stratification Model"""

# ╔═╡ 52cc571c-7d29-4167-834c-6672bbef6edf
md"""**Define Typed Stratification Model**"""

# ╔═╡ 1561268d-bec2-4b10-8f7f-996226a38e43
begin
	mdl_strat = LabelledPetriNet([:Q,:NQ],
    	:quarantine => ((:NQ)=>(:Q)),
    	:unquarantine => ((:Q)=>(:NQ)));
end;

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
	strat_ss = StrataSpec(mdl_strat_typed, [[:disease], [:disease,:infect]])
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
md"""**SIRD-Mask Model**"""

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

# ╔═╡ c33580a9-1d22-42a5-b527-018c2448b209
begin
	SIRD_M_ss = StrataSpec(SIRD_typed, [[:strata],[:strata],[:strata],[]])
	Mask_ss = StrataSpec(Mask_typed, [[:disease,:infect], [:disease,:infect]])
	mdl_SIRD_M, obsSIRD_M = stratify(SIRD_M_ss, Mask_ss, types′);
end;

# ╔═╡ a9043c9e-ec1b-41e8-a617-5a1a4d611415
AlgebraicPetri.Graph(mdl_SIRD_M)

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
begin
	true_mdl = mdl_SIRD_Q
	true_obs = obsSIRD_Q
	true_p = [5.0e-4, 5.0e-4, 5.0e-4, 5.0e-4, 5.0e-4, 5.0e-4,
            1.0e-2, 1.0e-5, 1.0e-3, 1.0e-5, 5.0e-5]
	u0 = [0.0,0.0,0.0,0.0,999.0,1.0,0.0,0.0] 
	tspan = (0.0,250.0)
	sample_data, sample_times, prob_real, true_sol, noiseless_data, data_labels = 		generateData(true_mdl, true_p, u0, tspan, 50, true_obs)
end;

# ╔═╡ 13f557c7-9002-4a2d-8be5-fccef5aa9dc6
md"""### Calibrate"""

# ╔═╡ 3b8954ae-2959-4ced-bdd2-1381121eb6c9
begin
	states_to_count = [:S,:I,:R,:D]
	p_init = repeat([1.0e-6], nt(true_mdl))
	p_est, sol_est, loss = calibrate(true_mdl, true_obs, states_to_count, u0, p_init, 	sample_data, sample_times, data_labels)
end;

# ╔═╡ 2276bfa9-c6bd-4c54-a6e1-7a88fd7a7762
md"""### Results"""

# ╔═╡ 79633d27-71b4-4261-81dc-db609133e2cd
md"""**Estimated Parameters**"""

# ╔═╡ ef69c329-9e2d-444a-8fa5-447f15b419b1
hcat(true_mdl[:tname],p_est)

# ╔═╡ 60408568-5762-43aa-9fef-b10a8c903422
md"""**Estimated Observations**"""

# ╔═╡ 8f1bfc19-0c0d-4ac9-a16a-fb56555f6707
plot_obs_w_ests(sample_times, sample_data, sol_est, true_obs)

# ╔═╡ Cell order:
# ╠═915f4df0-bc3a-11ec-2560-c9772e143679
# ╠═60d7decb-a7fa-494e-93ca-26d8b957cfcf
# ╠═aad0407a-4979-4ea4-94d3-38f1bb6d3de1
# ╟─287a44bb-e4d9-4a8b-b383-13d9c01e6e1e
# ╟─4330976d-ceeb-4a75-97d0-a8375ede795b
# ╠═80e375d8-60fc-4857-b430-9973117c5d29
# ╠═e2c11399-7878-4602-a276-e190857b3fa6
# ╟─2c24347b-9c55-4e14-876a-ea8e2867eaa2
# ╠═f0767520-04f6-4e97-a889-cf5e45d70b4c
# ╟─fbe5adf1-1d6b-40a6-b36e-49e3fc76057e
# ╠═4c3d24e2-254e-427e-9eb5-2344c1e8406d
# ╟─12a51952-34ab-47a3-a318-c4f58e3f3c12
# ╟─52cc571c-7d29-4167-834c-6672bbef6edf
# ╠═1561268d-bec2-4b10-8f7f-996226a38e43
# ╠═5b66805c-42ef-4611-a31b-962ba3f549ed
# ╠═b598e5b8-38e2-4d7d-afc9-2b0e7b1d2be8
# ╟─ec24996f-12f8-4d82-b210-8b23529fa73b
# ╠═be5b29e0-acb5-4bdd-951e-26fb1150a827
# ╟─1a35abc0-66cf-41ae-b313-7dd248cef669
# ╠═2f2d6daf-c0b8-4501-86a5-65a844c9359e
# ╠═dad0382a-3656-4bc9-acf6-b914916b15ac
# ╟─1c42f5e2-4e5c-4616-ab2d-c4d3ac877b9b
# ╠═cbcc0224-dae5-49c3-81eb-b1d37abbd526
# ╟─6b2b41bd-ea99-4917-8eec-d6963243f35a
# ╠═1a5c2273-a1d1-44a4-a312-db8e9cd3c7de
# ╠═44328cee-fbf9-42db-afeb-7fc48815952e
# ╠═84fec776-b42c-4066-ad45-624b1eb93460
# ╠═d482440e-2cce-4f3d-b4cb-208861a71187
# ╠═f7d18d80-f2a6-474c-8115-f323e09637d8
# ╠═c4f93cf2-9b30-44b2-9b7a-beee37686f37
# ╠═3c922e81-61af-4789-b2f4-e96dd981810f
# ╠═cc4d55cc-3ad1-47ed-9755-2ef3c616bf57
# ╠═0de945ad-093d-4975-8e48-a0748e2414b1
# ╠═385a2551-4733-4fb0-9d5d-9763f6757578
# ╠═452ade8f-6140-4531-9919-51e899b15f33
# ╠═c33580a9-1d22-42a5-b527-018c2448b209
# ╠═a9043c9e-ec1b-41e8-a617-5a1a4d611415
# ╠═d9fd32c2-ed42-4af6-90a7-d144296e222d
# ╠═4c7d9df5-8f92-4fe1-a2e3-82f607d28ccb
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
