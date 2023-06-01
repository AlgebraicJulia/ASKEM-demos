### A Pluto.jl notebook ###
# v0.19.16

using Markdown
using InteractiveUtils

# ╔═╡ e7806d7a-70fa-11ed-0d74-b9ab3ac85fa1
begin
	using Pkg
	Pkg.activate("..")
	using Revise, ASKEM
	using DifferentialEquations
	using Plots

	using Catlab.Theories
	using Catlab.Graphics
	using Catlab.CategoricalAlgebra
	using Catlab.Programs.RelationalPrograms

	using AlgebraicPetri
	using AlgebraicPetri.Epidemiology

	using ModelingToolkit
	using PlutoUI
	using HypertextLiteral
end

# ╔═╡ 40d6e607-0223-4b76-a262-034f3b1418bf


# ╔═╡ Cell order:
# ╠═e7806d7a-70fa-11ed-0d74-b9ab3ac85fa1
# ╠═40d6e607-0223-4b76-a262-034f3b1418bf
