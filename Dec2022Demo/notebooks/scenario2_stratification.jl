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

	using ASKEM.Dec2022Demo: formSIRD, formInfType, augLabelledPetriNet, sirdAugStates, typeSIRD, makeMultiAge, typeAge, typed_stratify, formVax, vaxAugStates, typeVax, writeMdlStrat, draw
	using ASKEM.Upstream: presentationToLabelledPetriNet, deserialize_wiringdiagram

end

# ╔═╡ 1df5e4a8-e761-4767-aed4-866534398922
md"""## Scenario 2: Stratifying by Age and Vaccination"""

# ╔═╡ 98ba5842-8ff9-4bce-a0a1-04c228b9ff9c
md"""Here we take a base SIRD disease model and stratify it by:
- 7 age groups and 
- vaccination status
"""

# ╔═╡ 943d7c37-ade5-4a4d-ae0b-e2bf390fcf6f
md"""### SIRD Model"""

# ╔═╡ 85ddcf05-c9c5-4436-8511-3713dbbe7694
begin
	SIRD = read_json_acset(LabelledPetriNet,"../SIRD.json")
	AlgebraicPetri.Graph(SIRD)
end

# ╔═╡ c7535760-7c2d-4552-a39a-7b51af04a92f
md"""### Age Model"""

# ╔═╡ 61a6b383-b9e2-4acf-8bb7-0fa926a0ec12
md"""### Vaccination Model"""

# ╔═╡ b323ac71-837a-4da1-ac81-05ac1ca9d600
begin
	StratificationWorkflow_lpn = read_json_acset(LabelledPetriNet,"../s2_strat_wf_present.json")
	AlgebraicPetri.Graph(StratificationWorkflow_lpn)
end

# ╔═╡ 6762c0ea-a453-451c-80c1-566a1389a45c
begin
	stratify_sird_age_vax = deserialize_wiringdiagram("../s2_strat_sird_age_vax.json")
	draw(stratify_sird_age_vax)
end

# ╔═╡ bb368013-7b99-4b2d-84e1-d22466983a74
begin 
	stratify_sird_hom_expr = to_hom_expr(FreeBiproductCategory, stratify_sird_age_vax)
	stratify_sird_jfunc = Catlab.Programs.GenerateJuliaPrograms.compile_expr(stratify_sird_hom_expr)
end;

# ╔═╡ 80bb3bb1-03e9-450c-ba96-18ea19a29a6d
eval(stratify_sird_jfunc)(7,"sird_age7_vax.json")

# ╔═╡ b602a85d-6792-4d93-b050-e901de024f5d
AlgebraicPetri.Graph(dom(SIRD_AGE_Vax))

# ╔═╡ 505cd6b7-f26c-42e0-812b-6962255d3648
md"""### SVIIvR Disease Model"""

# ╔═╡ 7f5c4e33-1c8b-41a1-ba73-52d3643fa008


# ╔═╡ Cell order:
# ╠═9aec393e-a083-45f5-ad73-e6bef22bb056
# ╠═1df5e4a8-e761-4767-aed4-866534398922
# ╠═98ba5842-8ff9-4bce-a0a1-04c228b9ff9c
# ╠═943d7c37-ade5-4a4d-ae0b-e2bf390fcf6f
# ╠═85ddcf05-c9c5-4436-8511-3713dbbe7694
# ╠═c7535760-7c2d-4552-a39a-7b51af04a92f
# ╠═61a6b383-b9e2-4acf-8bb7-0fa926a0ec12
# ╠═b323ac71-837a-4da1-ac81-05ac1ca9d600
# ╠═6762c0ea-a453-451c-80c1-566a1389a45c
# ╠═bb368013-7b99-4b2d-84e1-d22466983a74
# ╠═80bb3bb1-03e9-450c-ba96-18ea19a29a6d
# ╠═b602a85d-6792-4d93-b050-e901de024f5d
# ╠═505cd6b7-f26c-42e0-812b-6962255d3648
# ╠═7f5c4e33-1c8b-41a1-ba73-52d3643fa008
