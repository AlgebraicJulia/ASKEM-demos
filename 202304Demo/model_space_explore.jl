using Catlab, AlgebraicPetri, ModelingToolkit, DifferentialEquations, UnPack, SciMLBase,
      Distributions, Symbolics, DiffEqBase, Plots, EasyModelAnalysis
using DifferentialEquations.EnsembleAnalysis
using Catlab.CategoricalAlgebra
# using OpenAIReplMode
using CSV, DataFrames

using AlgebraicPetri.SubACSets

dir = "/Users/stokes.p/projects/AlgebraicJulia/ASKEM-demos/202304Demo/mdls"
fns = filter(endswith("json"), readdir(dir; join = true))
sch_fn = "/Users/stokes.p/projects/AlgebraicJulia/py-acsets/src/acsets/schemas/catlab/PropertyLabelledReactionNet.json"

SchPLRN = read_json_acset_schema(sch_fn)

@acset_type LPRN3(SchPLRN) <: AbstractLabelledReactionNet
miras = []
for fn in fns
    push!(miras, read_json_acset(LPRN3{Symbol, Symbol, Any, Any}, fn))
end
m = miras[1]


#*******

function strip_names(p::AbstractLabelledReactionNet)
    map(p, Name = name -> nothing)
  end
  

tmp1 = LabelledReactionNet{Any,Any}()
copy_parts!(tmp1,miras[1])
m1 = LabelledPetriNet(tmp1)

tmp2 = LabelledReactionNet{Any,Any}()
copy_parts!(tmp2,miras[2])
m2 = LabelledPetriNet(tmp2)

tmp3 = LabelledReactionNet{Any,Any}()
copy_parts!(tmp3,miras[3])
m3 = LabelledPetriNet(tmp3)

m1b = PetriNet()
copy_parts!(m1b,miras[1])

m3b = PetriNet()
copy_parts!(m3b,miras[3])


test = mca(m1,m2)
boom = []
for ii in test 
    boom = mca(test,m3)
end

#****



