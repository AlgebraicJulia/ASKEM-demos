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

function wf_load(f)
    return read_json_acset(LabelledReactionNet{Float64,Float64},f)
end

function wf_get_dim(m)
    return nt(m)
end

function wf_mcopy(m)
    return pair(m,m)
end

function wf_get_rates(m)
    return m[:rate]
end

function wf_rand(s,num_s)
    rng = Random.seed!(s)
    v = rand(rng,num_s)
    return v
end

function wf_rate_add(m,v)
    m[:rate] = m[:rate] + v 
    return m
end

# @present Workflow(FreeSymmetricMonoidalCategory) begin
@present Workflow(FreeBiproductCategory) begin
    (File,LRN,Dim,Seed,RVect,TSpan,CExpr)::Ob 
    wf_load::Hom(File,LRN)
    # wf_mcopy::Hom(LRN,LRN⊗LRN)
    wf_get_dim::Hom(LRN,Dim)
    wf_get_rates::Hom(LRN,RVect)
    wf_rand::Hom(Seed⊗Dim,RVect)
    wf_rate_add::Hom(LRN⊗RVect,LRN)
    MTKCompile::Hom(LRN⊗TSpan,CExpr)
end

load_perturb_sim = @program Workflow (f::File,s::Seed,ts::TSpan) begin # 
    m = wf_load(f)
    # m1, m2 = wf_mcopy(m)
    n_param = wf_get_dim(m)
    v = wf_rand(s,n_param)
    # v = wf_get_rates(m)
    m_perturb = wf_rate_add(m,v)
    sim_expr = MTKCompile(m_perturb,ts)
    return sim_expr # m_perturb # 
end

make_depvar(p,t) = :($p($t))
# function MTKCompile(bn::Union{AbstractLabelledBilayerNetwork, AbstractBilayerNetwork})
function MTKCompile(lrxn::AbstractLabelledReactionNet,tspan)    
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

    # sys = ODESystem(eq, t, name=:PetriNet)
    # prob = ODEProblem(sys, c, tspan, r)
    # sol = solve(prob,Tsit5())
  
    # prob = ODEProblem(MakeReactionSystem(model), u0, tspan, p)
    # sol = solve(prob, Tsit5(), tstops=sample_times)

    return quote
      $varstmt
      $paramstmt
      $diffstmt
      $(ϕs...)
      $eqnstmt
      return $varnameexpr, $parnameexpr, solve(ODEProblem(ODESystem(eqs, t, name=:PetriNet), $c, $tspan, $r))
    end
  
end
  
# serialize program wiring diagram
write_json_acset(load_perturb_sim.diagram, "diagram.json")

# visualize simulation plan
draw(load_perturb_sim)

# generate Julia program that executes simulation plan
wf_hom_expr = to_hom_expr(FreeBiproductCategory,load_perturb_sim)
wf_jfunc = Catlab.Programs.GenerateJuliaPrograms.compile(wf_hom_expr)

# expected output
#=  === wf_jfunc == 
function = (x1, x2, x3;) -> begin
    begin
        v1 = (Main).wf_load(x1)
        v2 = (Main).wf_get_dim(v1)
        v3 = (Main).wf_rand(x2, v2)
        v4 = (Main).wf_rate_add(v1, v3)
        v5 = (Main).MTKCompile(v4, x3)
        return v5
    end
end
=#

# apply plan to lrxnet from MIRA integration demo
wf_script = wf_jfunc(joinpath(@__DIR__, "..", "Oct2022Demo", "lrxnet_Mira_TC_est.json"),1234,(0,50))

# Run that simulation!
# wf_vars, wf_params, wf_ode_sol = eval(wf_script);