using ModelingToolkit
using AlgebraicPetri
using AlgebraicPetri.Epidemiology
using AlgebraicPetri.BilayerNetworks

using Catlab
using Catlab.CategoricalAlgebra
using Catlab.Graphics
import Catlab.CategoricalAlgebra: migrate!
using Catlab.WiringDiagrams
using Catlab.Programs.RelationalPrograms
display_uwd(ex) = to_graphviz(ex, box_labels=:name, junction_labels=:variable, edge_attrs=Dict(:len=>".75"));

sir = @relation (s,i,r) begin
    infection(s,i)
    recovery(i,r)
end
display_uwd(sir)

psir = apex(oapply_epi(sir))
psir
Graph(psir)

bnsir = LabelledBilayerNetwork()
migrate!(bnsir, psir)
bnsir
to_graphviz(bnsir)

make_depvar(p,t) = :($p($t))

@assert make_depvar(:a, :t) == :(a(t))
@assert make_depvar(:y, :x) == :(y(x))

function compile(bn::Union{AbstractLabelledBilayerNetwork, AbstractBilayerNetwork})
  varstmt = :(@variables t)
  @show varnames = bnsir[:variable]
  append!(varstmt.args, make_depvar.(bnsir[:variable], :t))

  paramstmt = :(@parameters)
  params = bnsir[:parameter]
  append!(paramstmt.args, bnsir[:parameter])

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

  return quote
    $varstmt
    $paramstmt
    $diffstmt
    $(ϕs...)
    $eqnstmt
    return $varnameexpr, $parnameexpr, ODESystem(eqs, t, name=:PetriNet)
  end

end

fex = compile(bnsir)

quote
  #= /Users/fairbanksj/github/AlgebraicJulia/ASKEM-demos/MTK/petri.jl:39 =# 
  @variables t S(t) I(t) R(t)
  @parameters inf rec
  D = Differential(t)
  ϕ1 = inf * S * I
  ϕ2 = rec * I
  eqs = [D(S) ~ +(-ϕ1), D(I) ~ ϕ1 + ϕ1 + -ϕ1 + -ϕ2, D(R) ~ +ϕ2]
  return ODESystem(eqs, t, name = :PetriNet)
end
println(fex)

(params...)

begin
  @variables t S(t) I(t) R(t)
  @parameters inf rec
  D = Differential(t)
  ϕ1 = inf * S * I
  ϕ2 = rec * I
  eqs = [D(S) ~ +(-ϕ1), D(I) ~ ϕ1 + ϕ1 + -ϕ1 + -ϕ2, D(R) ~ +ϕ2]
  ((S, I, R), (inf, rec), ODESystem(eqs, t, name = :PetriNet))
end