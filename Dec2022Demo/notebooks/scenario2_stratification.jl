### A Pluto.jl notebook ###
# v0.19.12

using Markdown
using InteractiveUtils

# ╔═╡ 9aec393e-a083-45f5-ad73-e6bef22bb056
begin 
	using Pkg#, Revise
	Pkg.activate(Base.current_project())
	
	using Catlab, Catlab.Theories
	using Catlab.CategoricalAlgebra
	using Catlab.Graphics
	using Catlab.Graphics: Graphviz
	import Catlab.CategoricalAlgebra: migrate!
	using Catlab.WiringDiagrams
	using Catlab.Programs
	using Catlab.Programs.RelationalPrograms
	
	using Catlab.Present
	using AlgebraicPetri
	using AlgebraicPetri: Graph

	using ASKEM.Dec2022Demo: formSIRD, formInfType, augLabelledPetriNet, sirdAugStates, typeSIRD, makeMultiAge, typeAge, typed_stratify, formVax, vaxAugStates, typeVax, writeMdlStrat, draw, loadSVIIvR, sviivrAugStates, typeSVIIvR
	using ASKEM.Upstream: presentationToLabelledPetriNet, deserialize_wiringdiagram
	using Base.Iterators
end

# ╔═╡ 1df5e4a8-e761-4767-aed4-866534398922
md"""## Scenario 2: Stratifying by Age and Vaccination"""

# ╔═╡ af1f77c5-6177-469e-a769-e0813aa61d42
md"""We approach the stratification two ways:
1. Stratifying a disease model (SIRD) by age (7 groups) and then by vaccination status
2. Stratifying a disease model already containing vax (SVIIvR) by age (7 groups)"""

# ╔═╡ 366eb04f-820e-4c09-a316-1d15405c981e
md"""The SVIIvR model used is read in from a pre-formed bilayer."""

# ╔═╡ d3e41be3-1cf7-469f-87a0-9296aa98e506
md"""We compute the stratifications here via wiring diagrams. To do so, we:
1. Define the presentation of objects and homomophisms to be used in the wiring diagrams.
2. Define the wiring diagrams of the individual stratification workflows.
3. Convert the wiring diagrams to expressions.
4. Compile the expressions into executable functions."""

# ╔═╡ 229e39cc-54f1-46b1-98a4-80d8536f2e42
md"""### SIRD Stratified by Age and Vax"""

# ╔═╡ 98ba5842-8ff9-4bce-a0a1-04c228b9ff9c
md"""Here we take a base SIRD disease model and stratify it by:
- 7 age groups and 
- vaccination status
"""

# ╔═╡ 943d7c37-ade5-4a4d-ae0b-e2bf390fcf6f
md"""#### SIRD Model"""

# ╔═╡ 85ddcf05-c9c5-4436-8511-3713dbbe7694
begin
	idirpath = joinpath(@__DIR__,"../outputs")
	SIRD = read_json_acset(LabelledPetriNet,joinpath(idirpath,"SIRD.json"))
	AlgebraicPetri.Graph(SIRD)
end

# ╔═╡ c7535760-7c2d-4552-a39a-7b51af04a92f
md"""#### Age Model"""

# ╔═╡ acba5303-3bee-41b6-b658-7c27bdc05691
md"""We display the age model for 3 groups, for ease of visualization, but stratify by 7 groups."""

# ╔═╡ 20fe503a-ddd2-4142-811a-0e6b81f13975
makeMultiAge(3) |> AlgebraicPetri.Graph

# ╔═╡ 61a6b383-b9e2-4acf-8bb7-0fa926a0ec12
md"""#### Vaccination Model"""

# ╔═╡ b2cbdf31-7de3-46e8-b8d1-bb0b8993c610
formVax() |> AlgebraicPetri.Graph

# ╔═╡ 77db4e0d-13cb-463e-8ac1-242f7ba8e2ed
md"""#### Presentation of Stratification Workflows"""

# ╔═╡ b323ac71-837a-4da1-ac81-05ac1ca9d600
begin
	StratificationWorkflow_lpn = read_json_acset(LabelledPetriNet,joinpath(idirpath,"s2_strat_wf_present.json"))
	AlgebraicPetri.Graph(StratificationWorkflow_lpn)
end

# ╔═╡ 621be410-ab57-4769-9fb0-f5f96e149ac5
md"""#### Wiring Diagram of Workflow Stratifying SIRD by Age and Vax"""

# ╔═╡ 6762c0ea-a453-451c-80c1-566a1389a45c
begin
	stratify_sird_age_vax = deserialize_wiringdiagram(joinpath(idirpath,"s2_strat_sird_age_vax.json"))
	draw(stratify_sird_age_vax)
end

# ╔═╡ c536d47b-65fb-47c5-8c05-c88b8ec38392
md"""#### Form Expression, Compile, and Run"""

# ╔═╡ bb368013-7b99-4b2d-84e1-d22466983a74
begin 
	stratify_sird_hom_expr = to_hom_expr(FreeBiproductCategory, stratify_sird_age_vax)
	stratify_sird_jfunc = Catlab.Programs.GenerateJuliaPrograms.compile_expr(stratify_sird_hom_expr)
end;

# ╔═╡ 80bb3bb1-03e9-450c-ba96-18ea19a29a6d
eval(stratify_sird_jfunc)(7,joinpath(idirpath,"sird_age7_vax.json"));

# ╔═╡ 833838e7-88ad-48eb-9263-5059374341fa
md"""#### Stratified Model"""

# ╔═╡ b602a85d-6792-4d93-b050-e901de024f5d
begin
	SIRD_Age_Vax_rt = read_json_acset(AlgebraicPetri.LabelledPetriNetUntyped{Vector},joinpath(idirpath,"sird_age7_vax.json"))
end;

# ╔═╡ d6daa705-50b7-4a85-afd0-95a6b6d49697
map(SIRD_Age_Vax_rt, Name=x->Symbol.(vcat(x...))) |> AlgebraicPetri.Graph

# ╔═╡ 20afdbac-f282-4833-8b99-3a69c4debe9d
md"""#### Stratified Parameter Names"""

# ╔═╡ 81c5c5ba-3106-4926-9b05-50f0332ebd61
md"""The state and transition names are now tuples of the component names."""

# ╔═╡ 1299ffa8-28b7-4606-abd3-7cac0de12327
#= ratenames = map(SIRD_Age_Vax_rt, Name=x->Symbol.(flatten(x))) |> x->x[:tname]
ratenames = map(SIRD_Age_Vax_rt, Name=x->Symbol.([x[1][1],x[1][2],x[2]])) |> x->x[:tname] =#

# ╔═╡ c8f7dd26-4474-4cbc-8236-5896cdf5b33a
ratenames = map(SIRD_Age_Vax_rt, Name=x->Symbol.(vcat(x...))) |> x->x[:tname]

# ╔═╡ c1790e5c-c89c-4da7-9fdf-d25fbc699104
for i in 1:10
	println(ratenames[i])
end

# ╔═╡ 4a852c91-3eaf-4828-8934-fb3c0036e611
md"""### SVIIvR Stratified by Age"""

# ╔═╡ c7fb7085-6619-409d-86eb-ca7a2bbd696f
md"""Alternatively, we can start with the SVIIvR model and stratify it the 7 age groups.
"""

# ╔═╡ 2b09a5f9-1f60-4f4f-95aa-663675cd16f5
md"""Here we read in the SVIIvR disease model from a bilayer.
"""

# ╔═╡ 2c5e9b57-8056-421a-be9d-200ffda73510
md"""#### SVIIvR Model"""

# ╔═╡ dd8a007a-3199-4e35-b3ec-733aea26ca4e
begin
	idirpath2 = joinpath(@__DIR__,"../../data")
	SVIIvR = loadSVIIvR(joinpath(idirpath2,"CHIME_SVIIvR_dynamics_BiLayer.json"))
	AlgebraicPetri.Graph(SVIIvR)
end

# ╔═╡ cf02a89a-d3a1-4a2c-9806-72d9ed052ba6
md"""#### Wiring Diagram of Workflow Stratifying SVIIvR by Age"""

# ╔═╡ f179d38d-d52a-4431-8fe0-31f6182375d5
begin
	stratify_sviivr_age = deserialize_wiringdiagram(joinpath(idirpath,"s2_strat_sviivr_age.json"))
	draw(stratify_sviivr_age)
end

# ╔═╡ 531e0ece-74b7-4322-a5da-9613b2b5b60e
md"""#### Form Expression, Compile, and Run"""

# ╔═╡ f14037d2-107a-4d9e-83e8-38b5864bafc7
begin 
	stratify_sviivr_hom_expr = to_hom_expr(FreeBiproductCategory, stratify_sviivr_age)
	stratify_sviivr_jfunc = Catlab.Programs.GenerateJuliaPrograms.compile_expr(stratify_sviivr_hom_expr)
end;

# ╔═╡ 991eb9df-0583-4d38-8eee-89dad7eef1e7
eval(stratify_sviivr_jfunc)(joinpath(idirpath2,"CHIME_SVIIvR_dynamics_BiLayer.json"),7,joinpath(idirpath,"sviivr_age7.json"));

# ╔═╡ e0e8e53b-402a-4f93-bd3e-b81271654f59
md"""#### Stratified Model"""

# ╔═╡ a7ff432d-e3ac-4f70-917d-6f7003824a72
SVIIvR_Age_rt = read_json_acset(AlgebraicPetri.LabelledPetriNetUntyped{Vector},joinpath(idirpath,"sviivr_age7.json"));

# ╔═╡ 7da72eca-b60f-40ee-91fd-da6102a8e8f0
map(SVIIvR_Age_rt, Name=x->Symbol.(x)) |> AlgebraicPetri.Graph

# ╔═╡ c796fc21-590e-4ada-bf69-3672f46b0eee
md"""#### Stratified Parameter Names"""

# ╔═╡ 01562d05-467c-423f-bbe1-5e84eb6d1695
begin
	sviivr_ratenames = map(SVIIvR_Age_rt, Name=x->Symbol.(x)) |> x->x[:tname]
	for i in 1:10
		println(sviivr_ratenames[i])
	end
end

# ╔═╡ Cell order:
# ╠═9aec393e-a083-45f5-ad73-e6bef22bb056
# ╟─1df5e4a8-e761-4767-aed4-866534398922
# ╟─af1f77c5-6177-469e-a769-e0813aa61d42
# ╟─366eb04f-820e-4c09-a316-1d15405c981e
# ╟─d3e41be3-1cf7-469f-87a0-9296aa98e506
# ╟─229e39cc-54f1-46b1-98a4-80d8536f2e42
# ╟─98ba5842-8ff9-4bce-a0a1-04c228b9ff9c
# ╟─943d7c37-ade5-4a4d-ae0b-e2bf390fcf6f
# ╠═85ddcf05-c9c5-4436-8511-3713dbbe7694
# ╟─c7535760-7c2d-4552-a39a-7b51af04a92f
# ╟─acba5303-3bee-41b6-b658-7c27bdc05691
# ╠═20fe503a-ddd2-4142-811a-0e6b81f13975
# ╟─61a6b383-b9e2-4acf-8bb7-0fa926a0ec12
# ╠═b2cbdf31-7de3-46e8-b8d1-bb0b8993c610
# ╟─77db4e0d-13cb-463e-8ac1-242f7ba8e2ed
# ╠═b323ac71-837a-4da1-ac81-05ac1ca9d600
# ╟─621be410-ab57-4769-9fb0-f5f96e149ac5
# ╠═6762c0ea-a453-451c-80c1-566a1389a45c
# ╟─c536d47b-65fb-47c5-8c05-c88b8ec38392
# ╠═bb368013-7b99-4b2d-84e1-d22466983a74
# ╠═80bb3bb1-03e9-450c-ba96-18ea19a29a6d
# ╟─833838e7-88ad-48eb-9263-5059374341fa
# ╠═b602a85d-6792-4d93-b050-e901de024f5d
# ╠═d6daa705-50b7-4a85-afd0-95a6b6d49697
# ╟─20afdbac-f282-4833-8b99-3a69c4debe9d
# ╟─81c5c5ba-3106-4926-9b05-50f0332ebd61
# ╟─1299ffa8-28b7-4606-abd3-7cac0de12327
# ╠═c8f7dd26-4474-4cbc-8236-5896cdf5b33a
# ╠═c1790e5c-c89c-4da7-9fdf-d25fbc699104
# ╟─4a852c91-3eaf-4828-8934-fb3c0036e611
# ╟─c7fb7085-6619-409d-86eb-ca7a2bbd696f
# ╟─2b09a5f9-1f60-4f4f-95aa-663675cd16f5
# ╟─2c5e9b57-8056-421a-be9d-200ffda73510
# ╠═dd8a007a-3199-4e35-b3ec-733aea26ca4e
# ╟─cf02a89a-d3a1-4a2c-9806-72d9ed052ba6
# ╠═f179d38d-d52a-4431-8fe0-31f6182375d5
# ╟─531e0ece-74b7-4322-a5da-9613b2b5b60e
# ╠═f14037d2-107a-4d9e-83e8-38b5864bafc7
# ╠═991eb9df-0583-4d38-8eee-89dad7eef1e7
# ╟─e0e8e53b-402a-4f93-bd3e-b81271654f59
# ╠═a7ff432d-e3ac-4f70-917d-6f7003824a72
# ╠═7da72eca-b60f-40ee-91fd-da6102a8e8f0
# ╟─c796fc21-590e-4ada-bf69-3672f46b0eee
# ╠═01562d05-467c-423f-bbe1-5e84eb6d1695
