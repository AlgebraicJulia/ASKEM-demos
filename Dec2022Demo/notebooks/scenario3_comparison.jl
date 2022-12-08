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
	
	using ASKEM
	using Test
	using ASKEM.Ontologies: infectious_ontology, strip_names
	using ASKEM.Dec2022Demo: formAugSIR, formAugSIRD, formAugSIRD2, formAugQuarantine, altTypeSIR, altTypeSIRD, altTypeSIRD2, altTypeQuarantine, formTarget, formModelList
	using ASKEM.Upstream: presentationToLabelledPetriNet, deserialize_wiringdiagram

end

# ╔═╡ 1df5e4a8-e761-4767-aed4-866534398922
md"""## Scenario 3: Structural Model Comparison"""

# ╔═╡ 98ba5842-8ff9-4bce-a0a1-04c228b9ff9c
md"""Here we search a list of models for a pair whose product equals a target model"""

# ╔═╡ d628f016-9181-49ec-8fd9-ab5c219b615f
md"""In this example the target model is the product of SIRD and Quarantine components"""

# ╔═╡ a33613cc-e896-4e7c-8d88-766a460d68c4
md"""### ComparisonWorkflow Presentaion"""

# ╔═╡ b323ac71-837a-4da1-ac81-05ac1ca9d600
begin
	ComparisonWorkflow_lpn = read_json_acset(LabelledPetriNet,"../s3_compar_wf_present.json")
	AlgebraicPetri.Graph(ComparisonWorkflow_lpn)
end

# ╔═╡ f96cc889-2ff1-48f6-bc54-cb2eb511adc2
md"""### Wiring Diagram of Program to Find SIRD and Q Components of Product Model in List"""

# ╔═╡ 6762c0ea-a453-451c-80c1-566a1389a45c
begin
	find_sird_q_components = deserialize_wiringdiagram("../s3_find_sird_q.json")
	draw(find_sird_q_components)
end

# ╔═╡ c6c43847-6e3b-46d6-ac27-8e72c030f312
md"""Most of this program is the generation of the models. The search is within decompose."""

# ╔═╡ 2e319d94-bbc4-43ed-9f09-e30c0edde13f
md"""### Form Wiring Diagram into Expression and Compile"""

# ╔═╡ bb368013-7b99-4b2d-84e1-d22466983a74
begin 
	find_comps_hom_expr = to_hom_expr(FreeBiproductCategory, find_sird_q_components)
	find_comps_jfunc = Catlab.Programs.GenerateJuliaPrograms.compile_expr(find_comps_hom_expr)
end;

# ╔═╡ 598abb29-8e94-47cf-980e-4f880e47cd32
md"""### Run Compiled Expression"""

# ╔═╡ 80bb3bb1-03e9-450c-ba96-18ea19a29a6d
res = eval(find_comps_jfunc)()

# ╔═╡ 5b4e406e-7bed-49ca-a224-a73be0335618
md"""### Confirm Correct Components Identified """

# ╔═╡ d0ae3242-f36c-4383-8045-d438e8fa1b1b
md"""### Plot Identified Components """

# ╔═╡ b602a85d-6792-4d93-b050-e901de024f5d
AlgebraicPetri.Graph(dom(SIRD_AGE_Vax))

# ╔═╡ 505cd6b7-f26c-42e0-812b-6962255d3648
md"""### Maximal Common Subgraph"""

# ╔═╡ 7f5c4e33-1c8b-41a1-ba73-52d3643fa008


# ╔═╡ Cell order:
# ╠═9aec393e-a083-45f5-ad73-e6bef22bb056
# ╠═1df5e4a8-e761-4767-aed4-866534398922
# ╠═98ba5842-8ff9-4bce-a0a1-04c228b9ff9c
# ╠═d628f016-9181-49ec-8fd9-ab5c219b615f
# ╠═a33613cc-e896-4e7c-8d88-766a460d68c4
# ╠═b323ac71-837a-4da1-ac81-05ac1ca9d600
# ╠═f96cc889-2ff1-48f6-bc54-cb2eb511adc2
# ╠═6762c0ea-a453-451c-80c1-566a1389a45c
# ╠═c6c43847-6e3b-46d6-ac27-8e72c030f312
# ╠═2e319d94-bbc4-43ed-9f09-e30c0edde13f
# ╠═bb368013-7b99-4b2d-84e1-d22466983a74
# ╠═598abb29-8e94-47cf-980e-4f880e47cd32
# ╠═80bb3bb1-03e9-450c-ba96-18ea19a29a6d
# ╠═5b4e406e-7bed-49ca-a224-a73be0335618
# ╠═d0ae3242-f36c-4383-8045-d438e8fa1b1b
# ╠═b602a85d-6792-4d93-b050-e901de024f5d
# ╠═505cd6b7-f26c-42e0-812b-6962255d3648
# ╠═7f5c4e33-1c8b-41a1-ba73-52d3643fa008
