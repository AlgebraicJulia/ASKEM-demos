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

valueat(x::Number, u, t) = x
valueat(f::Function, u, t) = try f(u,t) catch e f(t) end
AlgebraicPetri.vectorfield(pn::AbstractPetriNet) = begin
    tm = TransitionMatrices(pn)
    dt = tm.output - tm.input
    f(du,u,p,t) = begin
      rates = zeros(eltype(du),nt(pn))
      # u_m = [u[sname(pn, i)] for i in 1:ns(pn)]
      # p_m = [p[tname(pn, i)] for i in 1:nt(pn)]
      u_m = u
      p_m = p
      for i in 1:nt(pn)
        rates[i] = valueat(p_m[i],u,t) * prod(u_m[j] ^ tm.input[i,j] for j in 1:ns(pn))
      end
      for j in 1:ns(pn)
        # du[sname(pn, j)] = sum(rates[i] * dt[i,j] for i in 1:nt(pn); init = 0.0)
        du[j] = sum(rates[i] * dt[i,j] for i in 1:nt(pn); init = 0.0)
      end
      return du
    end
    return f
end

end
