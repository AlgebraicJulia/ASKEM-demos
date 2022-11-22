using ModelingToolkit
using AlgebraicPetri
using AlgebraicPetri.Epidemiology
using AlgebraicPetri.BilayerNetworks

using Catlab, Catlab.Theories
using Catlab.CategoricalAlgebra
using Catlab.Graphics
using Catlab.Graphics: Graphviz
import Catlab.CategoricalAlgebra: migrate!
using Catlab.WiringDiagrams
using Catlab.Programs
using Catlab.Programs.RelationalPrograms
import Catlab.WiringDiagrams.DirectedWiringDiagrams: WiringDiagramACSet
import Catlab.CategoricalAlgebra.CSets: parse_json_acset

# using JSON

using Random
using DifferentialEquations

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

function MTKLoadLRN(f)
    return read_json_acset(LabelledReactionNet{Float64,Float64},f)
end

make_depvar(p,t) = :($p($t))
function MTKFormODEProb(lrxn::AbstractLabelledReactionNet,tspan)    
    r = lrxn[:rate]
    c = lrxn[:concentration]

    lpn = LabelledPetriNet(lrxn);
    bn = LabelledBilayerNetwork();
    migrate!(bn,lpn);
    
    varstmt = :(@variables t)
    @show varnames = bn[:variable]
    append!(varstmt.args, make_depvar.(bn[:variable], :t))
    
    paramstmt = :(@parameters)
    params = bn[:parameter]
    append!(paramstmt.args, bn[:parameter])
    
    diffstmt = :(D = Differential(t))
  
    ϕs = map(parts(bn, :Box)) do b
      vars = map(incident(bn, b,:call)) do i
        j = bn[i, :arg]
        return bn[j, :variable]
      end
      p = :(*($(bn[b, :parameter])))
      append!(p.args, vars)
      return :($(Symbol("ϕ$b")) = $p)
    end
  
    
    infs = map(parts(bn, :Qout)) do tv
      vars = map(incident(bn, tv, :infusion)) do wa
        j = bn[wa, :influx]
        return Symbol("ϕ$j")
      end
      p = :(+())
      append!(p.args, vars)
  
      # same for the outfluxes
      vars = map(incident(bn, tv, :effusion)) do wn
        j = bn[wn, :efflux]
        return :(- $(Symbol("ϕ$j")))
      end
      append!(p.args, vars)
      return p
    end
  
    zparts = zip(bn[:tanvar], infs)
  
    eqns = [:(D($tanvar) ~ $rhs) for (tanvar, rhs) in zparts]
    eq = :([])
    append!(eq.args, eqns)
    eqnstmt = :(eqs = $eq)
  
    varnameexpr = Expr(:tuple, varnames...)
    parnameexpr = Expr(:tuple, params...)

    eval(varstmt)
    eval(paramstmt)
    eval(diffstmt)
    map(eval,ϕs)
    eval(eqnstmt)
    sys = ODESystem(eqs, t, name=:PetriNet)
    prob = ODEProblem(sys, c, tspan, r)
    # sol = solve(prob,Tsit5())
  
    # return quote
    #   $varstmt
    #   $paramstmt
    #   $diffstmt
    #   $(ϕs...)
    #   $eqnstmt
    #   return $varnameexpr, $parnameexpr, ODEProblem(ODESystem(eqs, t, name=:PetriNet), $c, $tspan, $r)
    # end
    return prob
end

function MTKSolveODE(prob)
    return solve(prob,Tsit5())
end
  
# Form Workflow presentation of FreeBiproductCategory
@present Workflow(FreeBiproductCategory) begin
    (File,LRN,TSpan,ODEProb,ODESol)::Ob 
    MTKLoadLRN::Hom(File,LRN)
    MTKFormODEProb::Hom(LRN⊗TSpan,ODEProb)
    MTKSolveODE::Hom(ODEProb,ODESol)
end

# Form wiring diagram of load_form_sim Workflow
load_form_sim = @program Workflow (f::File,ts::TSpan) begin # 
    lrn = MTKLoadLRN(f)
    ode_prob = MTKFormODEProb(lrn,ts)
    ode_sol = MTKSolveODE(ode_prob)
    return ode_sol 
end

# Serialize program wiring diagram
# write_json_graph(load_form_sim,"diagram_load_form_sim.json") 
write_json_acset(load_form_sim.diagram, "diagram_load_form_sim.json")

# Example of reading in wiring diagram from JSON
# load_form_sim_roundtrip = read_json_graph(Symbol, Symbol, Nothing, "diagram_load_form_sim.json")
function parse_json_acset(::Type{T}, input::AbstractDict) where T <: ACSet # Catlab.CategoricalAlgebra.CSets.
  out = T()
  #=for (k,v) ∈ input
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
  end=#
  for (k,v) ∈ input
    add_parts!(out, Symbol(k), length(v))
  end
  for l ∈ values(input)
    for (i, j) ∈ enumerate(l)
      for k ∈ keys(j) # (k, v) = j # 
        v = j[k]
        vtype = eltype(out[Symbol(k)])
        # if ((Symbol(k)==:box_type) | (Symbol(k)==:outer_in_port_type) | (Symbol(k)==:inner_in_port_type)) & (v isa String)
        if v isa String
            v = Meta.parse(v)
            if v isa Expr
                v = eval(v)
            end
        end
        if !(v isa vtype)
          v = vtype(v)
        end
        set_subpart!(out, i, Symbol(k), v)
      end
    end
  end
  out
end

# Read in wiring diagram acset from file
rt_wd_acset = read_json_acset(WiringDiagramACSet{Any,Any,Any,DataType},"diagram_load_form_sim.json")
# rt_wd_acset = parse_json_acset(WiringDiagramACSet{Any,Any,Any,DataType},JSON.parsefile("diagram_load_form_sim.json"))

# Check equality of read-in wd-acset to original
rt_wd_acset == load_form_sim.diagram

# Form roundtrip wiring diagram from read-in wd acset
load_form_sim_roundtrip = WiringDiagram{ThBiproductCategory,Any,Any,Any}(rt_wd_acset,nothing)

# Check equality of roundtrip wiring diagram to original
load_form_sim == load_form_sim_roundtrip

# Visualize simulation plan
draw(load_form_sim_roundtrip)

# Generate Julia function that executes simulation plan
mtk_hom_expr = to_hom_expr(FreeBiproductCategory,load_form_sim_roundtrip)
mtk_jfunc = Catlab.Programs.GenerateJuliaPrograms.compile(mtk_hom_expr)

# Expected output
#=  === mtk_jfunc == 
function = (x1, x2;) -> begin
    begin
        v1 = (Main).MTKLoadLRN(x1)
        v2 = (Main).MTKFormODEProb(v1, x2)
        v3 = (Main).MTKSolveODE(v2)
        return v3
    end
end
=#

# Apply generated function to lrxnet from MIRA integration demo
mtk_ode_sol = mtk_jfunc(joinpath(@__DIR__, "..", "Oct2022Demo", "lrxnet_Mira_TC_est.json"),(0,50))
