module Upstream
export presentationToLabelledPetriNet, deserialize_wiringdiagram, draw

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
using Catlab.CategoricalAlgebra.CSets: parse_json_acset
function Catlab.CategoricalAlgebra.CSets.parse_json_acset(::Type{T}, input::AbstractDict) where T <: ACSet 
    out = T()
    for (k,v) ∈ input
      add_parts!(out, Symbol(k), length(v))
    end
    for l ∈ values(input)
      for (i, j) ∈ enumerate(l)
        for k ∈ keys(j) # (k,v) = j # 
          v = j[k]
          vtype = eltype(out[Symbol(k)])
          if !(v isa vtype)
            v = vtype(v)
          end
          set_subpart!(out, i, Symbol(k), v)
        end
      end
    end
    out
end

deserialize_wiringdiagram(filepath::String, value=nothing) = deserialize_wiringdiagram!(read_json_acset(WiringDiagramACSet{Any,Any,Any,Any}, filepath), isnothing(value) ? filepath : value)
deserialize_wiringdiagram!(dwd, value) = begin
  convsymbol(dwd, key) = begin
    dwd[key] .= Symbol.(dwd[key])
  end
  
  dwd[:box_type] .= Box{Symbol}
  convsymbol(dwd, :in_port_type)
  convsymbol(dwd, :out_port_type)
  convsymbol(dwd, :outer_in_port_type)
  convsymbol(dwd, :outer_out_port_type)
  dwd[:value] .= map(Symbol,dwd[:value])
  wd_acset2 = WiringDiagramACSet{Any,Any,Any,DataType}()
  copy_parts!(wd_acset2,dwd)
  return WiringDiagram{ThBiproductCategory, Any, Any, Any}(wd_acset2, value)
end


draw(d::WiringDiagram) = to_graphviz(d,
    orientation=LeftToRight,
    labels=true, label_attr=:xlabel,
    node_attrs=Graphviz.Attributes(
      :fontname => "Courier",
    ),
    edge_attrs=Graphviz.Attributes(
      :fontname => "Courier",
    )
)

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
