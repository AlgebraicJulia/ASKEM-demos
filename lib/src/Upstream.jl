module Upstream
export presentationToLabelledPetriNet

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

function presentationToLabelledPetriNet(present)
    lpn = LabelledPetriNet(map(Symbol,generators(present,:Ob)))
    state_idx = Dict()
    for ii in 1:ns(lpn)
        state_idx[sname(lpn,ii)] = ii
    end
    num_obs = length(generators(present,:Ob))
    num_homs = length(generators(present,:Hom))
    for curr_hom in generators(present,:Hom)
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
    return lpn
end

end