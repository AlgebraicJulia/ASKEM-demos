### A Pluto.jl notebook ###
# v0.19.13

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
	
	using ASKEM
	using Test
	using ASKEM.Ontologies: infectious_ontology, strip_names
	using ASKEM.Dec2022Demo: formAugSIR, formAugSIRD, formAugSIRD2, formAugQuarantine, altTypeSIR, altTypeSIRD, altTypeSIRD2, altTypeQuarantine, formTarget, formModelList
	using ASKEM.Upstream: presentationToLabelledPetriNet, deserialize_wiringdiagram

end

# ╔═╡ 7f5c4e33-1c8b-41a1-ba73-52d3643fa008
using Catlab.CategoricalAlgebra.CSets

# ╔═╡ 1df5e4a8-e761-4767-aed4-866534398922
md"""## Scenario 3: Structural Model Comparisons"""

# ╔═╡ a56072d9-c2a5-46ae-8139-8f7277a150f0
md"""We perform two types of structural model comparisons:
1. Search for components of a product model
2. Search for maximal common submodels
"""

# ╔═╡ b45e01a5-816d-48f6-9042-3ac3494e7925
md"""### Model Decomposition"""

# ╔═╡ 98ba5842-8ff9-4bce-a0a1-04c228b9ff9c
md"""Here we search a list of models for a pair whose product equals a target model"""

# ╔═╡ d628f016-9181-49ec-8fd9-ab5c219b615f
md"""In this example, the target model is the product of SIRD and Quarantine components"""

# ╔═╡ 39f07f7d-2ca0-462e-89d9-d7d11aece267
md"""The search list the component SIRD and Quarantine models in addition to SIR and  an alternative SIRD."""

# ╔═╡ a33613cc-e896-4e7c-8d88-766a460d68c4
md"""#### ComparisonWorkflow Presentaion"""

# ╔═╡ b323ac71-837a-4da1-ac81-05ac1ca9d600
begin
	idirpath_wd = joinpath(@__DIR__,"../outputs")
	ComparisonWorkflow_lpn = read_json_acset(LabelledPetriNet,joinpath(idirpath_wd,"s3_compar_wf_present.json"))
	AlgebraicPetri.Graph(ComparisonWorkflow_lpn)
end

# ╔═╡ f96cc889-2ff1-48f6-bc54-cb2eb511adc2
md"""#### Wiring Diagram of Program to Find SIRD and Quarantine Components of Product Model in a List"""

# ╔═╡ 6762c0ea-a453-451c-80c1-566a1389a45c
begin
	find_sird_q_components = deserialize_wiringdiagram(joinpath(idirpath_wd,"s3_find_sird_q.json"))
	draw(find_sird_q_components)
end

# ╔═╡ c6c43847-6e3b-46d6-ac27-8e72c030f312
md"""Most of this program is the generation of the models. The search is within decompose."""

# ╔═╡ 2e319d94-bbc4-43ed-9f09-e30c0edde13f
md"""#### Form Wiring Diagram into Expression and Compile"""

# ╔═╡ bb368013-7b99-4b2d-84e1-d22466983a74
begin 
	find_comps_hom_expr = to_hom_expr(FreeBiproductCategory, find_sird_q_components)
	find_comps_jfunc = Catlab.Programs.GenerateJuliaPrograms.compile_expr(find_comps_hom_expr)
end;

# ╔═╡ 598abb29-8e94-47cf-980e-4f880e47cd32
md"""#### Run Compiled Expression"""

# ╔═╡ 80bb3bb1-03e9-450c-ba96-18ea19a29a6d
res = eval(find_comps_jfunc)();

# ╔═╡ d0ae3242-f36c-4383-8045-d438e8fa1b1b
md"""#### Plot/Confirm Identified Components """

# ╔═╡ c8a76821-c342-43ea-9f68-29c507d35303
md"""The result comprises two pairs of models."""

# ╔═╡ d0476e7b-e8df-4b48-b51f-42235c65cab4
md"""The first pair is the SIRD and Quarantine."""

# ╔═╡ b602a85d-6792-4d93-b050-e901de024f5d
AlgebraicPetri.Graph(dom(res[1][1]))

# ╔═╡ d5db830a-4932-4c77-837e-26c4c917ee89
AlgebraicPetri.Graph(dom(res[1][2]))

# ╔═╡ 32be2bde-7b6a-49fb-b27a-779409ef3efd
md"""The second pair is Quarantine and the alternative SIRD."""

# ╔═╡ bee76011-267a-461b-91e4-2452cf940cee
dom(res[1][2]) == dom(res[2][1])

# ╔═╡ 6133a0c1-6b83-489f-b884-2fd558f66969
AlgebraicPetri.Graph(dom(res[2][2]))

# ╔═╡ d7efab3f-358a-4832-aee3-6f81ad9b5bff
is_isomorphic(strip_names(dom(res[1][1])), strip_names(dom(res[2][2])))

# ╔═╡ 505cd6b7-f26c-42e0-812b-6962255d3648
md"""### Maximal Common ACSet"""

# ╔═╡ 3a6ef240-1a22-497e-83b0-c7c0b01f58f5
md"""Suppose we are interested in two processes P₁ and P₂ which we suspect (perhaps thanks to our real-world inutitions) to be related, 
and suppose that we have modeled these processes with two models M₁ and M₂. Although it might be difficult to investigate the relationship between 
processes P₁ and P₂ (whose description may be mathematically opaque), it is significantly easier to investigate the relationship between Theories
mathematical models M₁ and M₂. Indeed, whenever M₁ and M₂ are represented as acsets, this investigation reduces to the problem of finding Maximum 
Common Submodel. This problem is a special case of the Maxium Common ACSet which is stated formally below."""

# ╔═╡ 0faae1c3-de39-466a-ad47-f07448e9a19c
md"""**Maximum Common Acset (MCA).**
Input: two Acsets a₁ and a₂.
Task : find all a with with |a| maximum possible such that there is a monic span of Acset a₁ ← a → a₂."""

# ╔═╡ 6048521c-ef8d-456a-8b56-c404ab61308f
md"""Finding a maximum common acset is easy in catlab: you can simply call mca(M₁, M₂). Let's demonstrate this. """

# ╔═╡ 17a5180d-784d-48cf-8e4c-eb00f42cb9cb
md"""#### MCA of SIR and SIS"""

# ╔═╡ c3811cb0-b1c1-470b-a05e-3fcabd3a5414
md"""For concreteness, let's work with Petri nets, specifically Petri nets which represent epidemiological models. 
One such net may be the **SIR model** shown below.  """

# ╔═╡ 791d1b45-0cf0-4fdd-870c-874ed13eb78b
begin
	idirpath_mdl = joinpath(@__DIR__,"../../data")
	sir = read_json_acset(LabelledPetriNet,joinpath(idirpath_mdl,"SIR.json"))
	AlgebraicPetri.Graph(sir)
end

# ╔═╡ aa56c6c4-2c33-48a8-8620-66f1a38271de
md"""Another such Petri net might be the **SIS model**. We've drawn it below."""

# ╔═╡ d4a49ab0-9a7c-4dfd-ae24-b31502a77ec2
begin
	sis    = read_json_acset(LabelledPetriNet,joinpath(idirpath_mdl,"SIS.json"))
	AlgebraicPetri.Graph(sis)
end

# ╔═╡ 6eaa9d5a-07e3-4c92-a73e-75e9d3a6fd85
md"""Now, if we are interested in determining all of the maximum common submodels of the SIR and SIS models, we can simply run"""

# ╔═╡ 6d975148-3259-427b-8146-8e030b07dc02
maxcommon = mca(sir, sis);

# ╔═╡ ebeba408-a32c-4534-a5b2-96103a42afde
md"""In this case, there are two maximum common submodels; these are:"""

# ╔═╡ f1ea45bc-8e43-4e74-be93-d705ba6e266e
AlgebraicPetri.Graph(maxcommon[2])

# ╔═╡ 3b20b669-26f6-4464-824e-7f9575e0682f
AlgebraicPetri.Graph(maxcommon[1])

# ╔═╡ Cell order:
# ╠═9aec393e-a083-45f5-ad73-e6bef22bb056
# ╟─1df5e4a8-e761-4767-aed4-866534398922
# ╟─a56072d9-c2a5-46ae-8139-8f7277a150f0
# ╟─b45e01a5-816d-48f6-9042-3ac3494e7925
# ╟─98ba5842-8ff9-4bce-a0a1-04c228b9ff9c
# ╟─d628f016-9181-49ec-8fd9-ab5c219b615f
# ╟─39f07f7d-2ca0-462e-89d9-d7d11aece267
# ╟─a33613cc-e896-4e7c-8d88-766a460d68c4
# ╠═b323ac71-837a-4da1-ac81-05ac1ca9d600
# ╟─f96cc889-2ff1-48f6-bc54-cb2eb511adc2
# ╠═6762c0ea-a453-451c-80c1-566a1389a45c
# ╟─c6c43847-6e3b-46d6-ac27-8e72c030f312
# ╟─2e319d94-bbc4-43ed-9f09-e30c0edde13f
# ╠═bb368013-7b99-4b2d-84e1-d22466983a74
# ╟─598abb29-8e94-47cf-980e-4f880e47cd32
# ╠═80bb3bb1-03e9-450c-ba96-18ea19a29a6d
# ╟─d0ae3242-f36c-4383-8045-d438e8fa1b1b
# ╟─c8a76821-c342-43ea-9f68-29c507d35303
# ╟─d0476e7b-e8df-4b48-b51f-42235c65cab4
# ╠═b602a85d-6792-4d93-b050-e901de024f5d
# ╠═d5db830a-4932-4c77-837e-26c4c917ee89
# ╟─32be2bde-7b6a-49fb-b27a-779409ef3efd
# ╠═bee76011-267a-461b-91e4-2452cf940cee
# ╠═6133a0c1-6b83-489f-b884-2fd558f66969
# ╠═d7efab3f-358a-4832-aee3-6f81ad9b5bff
# ╟─505cd6b7-f26c-42e0-812b-6962255d3648
# ╠═7f5c4e33-1c8b-41a1-ba73-52d3643fa008
# ╟─3a6ef240-1a22-497e-83b0-c7c0b01f58f5
# ╟─0faae1c3-de39-466a-ad47-f07448e9a19c
# ╟─6048521c-ef8d-456a-8b56-c404ab61308f
# ╟─17a5180d-784d-48cf-8e4c-eb00f42cb9cb
# ╟─c3811cb0-b1c1-470b-a05e-3fcabd3a5414
# ╠═791d1b45-0cf0-4fdd-870c-874ed13eb78b
# ╟─aa56c6c4-2c33-48a8-8620-66f1a38271de
# ╠═d4a49ab0-9a7c-4dfd-ae24-b31502a77ec2
# ╟─6eaa9d5a-07e3-4c92-a73e-75e9d3a6fd85
# ╠═6d975148-3259-427b-8146-8e030b07dc02
# ╟─ebeba408-a32c-4534-a5b2-96103a42afde
# ╠═f1ea45bc-8e43-4e74-be93-d705ba6e266e
# ╠═3b20b669-26f6-4464-824e-7f9575e0682f
