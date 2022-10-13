### A Pluto.jl notebook ###
# v0.19.13

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

# ╔═╡ 3eeb61e2-49d1-11ed-23db-c961288aede3
begin 
	using Pkg#, Revise
	Pkg.activate(Base.current_project())
	# Pkg.instantiate()
	using Catlab.CategoricalAlgebra
	using Catlab.Present, Catlab.Theories
	using AlgebraicPetri
	using AlgebraicPetri: Graph
	using Semagrams
	using Catlab.Graphics
end;

# ╔═╡ 69a52c7a-6877-4ec4-9bea-e215cde12c91
include("StratPetris.jl")

# ╔═╡ 0d9eff21-cf7b-4de1-befb-fd5437468acd
include("Oct2022Demo.jl")

# ╔═╡ d4f1d2f4-d69d-4e0c-b219-eecd08c44768
@bind mdl_strat_sema Semagram{StratPetri}("https://semagrams-builds.s3.amazonaws.com/aff1944/petri/main.js",
	"Petri",
	Dict{Symbol,Function}(
		:Name => s -> Symbol(s),
		:Rate => s -> s,
		:Concentration => s -> s,
		:TransitionTypeValue => s -> 
			if length(s) == 0
				:none
			else
				Symbol(s[1])
			end,
		:Stratification => s -> Symbol[Symbol(str) for str in s]
	),
)

# ╔═╡ 048a8b5f-2eea-4189-ab4b-9c2c757cf3d5
begin
	types′ = LabelledPetriNet([:Pop],
    	:infect=>((:Pop, :Pop)=>(:Pop, :Pop)),
    	:disease=>(:Pop=>:Pop),
    	:strata=>(:Pop=>:Pop))
	types = map(types′, Name=name->nothing)
end;

# ╔═╡ 6cc5c1bc-958c-4a70-8dda-99ea74a0f7e2
mdl_strat = begin
	p = LabelledPetriNet()
	copy_parts!(p, mdl_strat_sema)
	p
end;

# ╔═╡ a9eeaced-d364-461a-99fc-662025c32cc6
function gettype(acs, i)
	findfirst([:infect, :disease, :strata] .== subpart(acs, i, :transitiontype))
end

# ╔═╡ 77ddb8fb-8a6f-4528-b184-9729d163d197
inttypes = Int[gettype(mdl_strat_sema, i) for i in parts(mdl_strat_sema, :T)]

# ╔═╡ 91c6bb81-0ae4-4706-881d-d220613eff2d
	mdl_strat_typed = homomorphism(mdl_strat, types;
    	initial=(T=inttypes,), type_components=(Name=x->nothing,))

# ╔═╡ 73e95f2c-2766-4b8b-85f4-92f7c63ec8be
strat_ss = StrataSpec(mdl_strat_typed, subpart(mdl_strat_sema, :stratificationwith))

# ╔═╡ Cell order:
# ╠═3eeb61e2-49d1-11ed-23db-c961288aede3
# ╠═69a52c7a-6877-4ec4-9bea-e215cde12c91
# ╠═0d9eff21-cf7b-4de1-befb-fd5437468acd
# ╠═d4f1d2f4-d69d-4e0c-b219-eecd08c44768
# ╠═048a8b5f-2eea-4189-ab4b-9c2c757cf3d5
# ╠═6cc5c1bc-958c-4a70-8dda-99ea74a0f7e2
# ╠═a9eeaced-d364-461a-99fc-662025c32cc6
# ╠═77ddb8fb-8a6f-4528-b184-9729d163d197
# ╠═91c6bb81-0ae4-4706-881d-d220613eff2d
# ╠═73e95f2c-2766-4b8b-85f4-92f7c63ec8be
