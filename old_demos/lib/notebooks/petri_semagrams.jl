### A Pluto.jl notebook ###
# v0.19.16

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

# ╔═╡ 47d873d0-727b-11ed-0154-07418a359390
begin
	using Pkg
	Pkg.activate("..")
	using Revise, ASKEM
	using ASKEM.Semagrams.Petri
end

# ╔═╡ e57e6e95-727b-452d-8631-196738a1c04a
@bind p StratPetriSemagram()

# ╔═╡ 96209336-97c5-4dd2-a81a-29215ac51b12
typed_petri(p, infectious_ontology)

# ╔═╡ Cell order:
# ╠═47d873d0-727b-11ed-0154-07418a359390
# ╠═e57e6e95-727b-452d-8631-196738a1c04a
# ╠═96209336-97c5-4dd2-a81a-29215ac51b12
