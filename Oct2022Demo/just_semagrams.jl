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
end;

# ╔═╡ d4f1d2f4-d69d-4e0c-b219-eecd08c44768
@bind mdl_strat_sema Semagram{LabelledPetriNet}("https://semagrams-builds.s3.amazonaws.com/8282a34/petri/main.js",
	"Petri",
	Dict{Symbol,Function}(
		:Name => s -> Symbol(s),
		:Rate => s -> s,
		:Concentration => s -> s,
	),
)

# ╔═╡ 6abb6c42-cfd5-4aaa-8f89-1c53dde10e9a
mdl_strat_sema

# ╔═╡ Cell order:
# ╠═3eeb61e2-49d1-11ed-23db-c961288aede3
# ╠═d4f1d2f4-d69d-4e0c-b219-eecd08c44768
# ╠═6abb6c42-cfd5-4aaa-8f89-1c53dde10e9a
