using AlgebraicPetri

using Catlab, Catlab.Theories
using Catlab.CategoricalAlgebra
using Catlab.Graphics
using Catlab.Graphics: Graphviz
import Catlab.CategoricalAlgebra: migrate!
using Catlab.Graphs.BasicGraphs, Catlab.Graphs.BipartiteGraphs
using Catlab.WiringDiagrams
using Catlab.Programs
using Catlab.Programs.RelationalPrograms
import Catlab.WiringDiagrams.DirectedWiringDiagrams: WiringDiagramACSet


# Form Workflow presentation of FreeBiproductCategory
@present Workflow(FreeBiproductCategory) begin
    (File,LRN,TSpan,ODEProb,ODESol)::Ob 
    MTKLoadLRN::Hom(File,LRN)
    MTKFormODEProb::Hom(LRN⊗TSpan,ODEProb)
    MTKSolveODE::Hom(ODEProb,ODESol)
end

#**
# Construct Presentation as LabelledPetriNet
#**
lpn = LabelledPetriNet(map(Symbol,generators(Workflow,:Ob)))
state_idx = Dict()
for ii in 1:ns(lpn)
    state_idx[sname(lpn,ii)] = ii
end
num_obs = length(generators(Workflow,:Ob))
num_homs = length(generators(Workflow,:Hom))
for curr_hom in generators(Workflow,:Hom)
    i = add_transition!(lpn, tname=Symbol(curr_hom))
    num_args = length(args(dom(curr_hom)))
    if num_args==1
        add_inputs!(lpn,1,i,state_idx[Symbol(dom(curr_hom))])
    else
        add_inputs!(lpn,num_args,repeat([i],num_args),map(x->state_idx[Symbol(x)], args(dom(curr_hom))))
    end
    num_args = length(args(codom(curr_hom)))
    if num_args==1
        add_outputs!(lpn,1,i,state_idx[Symbol(codom(curr_hom))])
    else
        add_outputs!(lpn,num_args,repeat([i],num_args),map(x->state_idx[Symbol(x)], args(codom(curr_hom))))
    end
end

AlgebraicPetri.Graph(lpn)
write_json_acset(lpn,"presentation.json")

lpn_rt = read_json_acset(LabelledPetriNet,"presentation.json")
AlgebraicPetri.Graph(lpn_rt)

#**
# Construct Presentation as BipartiteGraph
#**
#=
num_obs = length(generators(Workflow,:Ob))
num_homs = length(generators(Workflow,:Hom))
bpg = BipartiteGraph()
for curr_ob in generators(Workflow,:Ob)
    add_vertex₁!(bpg) # ,label=Symbol(curr_ob))
end
for curr_hom in generators(Workflow,:Hom)
    add_vertex₂!(bpg) # ; name=Symbol(curr_hom))
end
for curr_hom in generators(Workflow,:Hom)
    for curr_ob in dom(curr_hom)

    end
    for curr_ob in codom(curr_hom)

    end
end
=#

#=
add_edge₁₂!(g, 1, 2)
add_edge₂₁!(g, 1, 1)
add_edges₁₂!(g, [2,2], [3,3])
add_edges₂₁!(g, [2,3], [1,1])
@test (ne₁₂(g), ne₂₁(g)) == (3,3)
@test (edges₁₂(g), edges₂₁(g)) == (1:3, 1:3)
@test ne(g) == (3,3)
@test edges(g) == (1:3, 1:3)
@test (src₁(g), tgt₂(g)) == ([1,2,2], [2,3,3])
@test (src₂(g), tgt₁(g)) == ([1,2,3], [1,1,1])

rem_edge₁₂!(g, 1, 2)
@test (src₁(g), tgt₂(g)) == ([2,2], [3,3])
rem_edge₂₁!(g, 1, 1)
@test (src₂(g), tgt₁(g)) == ([3,2], [1,1])
rem_vertex₁!(g, 1)
@test nv(g) == (1,3)
@test ne(g) == (2,0)
rem_vertex₂!(g, 3)
@test nv(g) == (1,2)
@test ne(g) == (0,0)
=#