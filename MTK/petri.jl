using Test
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
display_uwd(ex) = to_graphviz(ex, box_labels=:name, junction_labels=:variable, edge_attrs=Dict(:len => ".75"));

# Define SIR Model
sir = @relation (s, i, r) begin
  infection(s, i)
  recovery(i, r)
end
display_uwd(sir)

# Convert to Epidemiology petri net
psir = apex(oapply_epi(sir))
psir
Graph(psir)

# Create empty bilayer network
bnsir = LabelledBilayerNetwork()
# migrate petri model to bilayer network
migrate!(bnsir, psir)
bnsir
to_graphviz(bnsir)

make_depvar(p, t) = :($p($t))

@assert make_depvar(:a, :t) == :(a(t))
@assert make_depvar(:y, :x) == :(y(x))

function ModelingToolkit.ODESystem(bn::Union{AbstractLabelledBilayerNetwork,AbstractBilayerNetwork}; name = :PetriNet)
  t = (@variables t)[1]
  D = Differential(t)
  symbolic_vars = map(bnsir[:variable]) do v
      (@variables $v(t))[1]
  end
  symbolic_params = map(bnsir[:parameter]) do p
      (@parameters $p)[1]
  end

  ϕs = map(parts(bn, :Box)) do b
    p = symbolic_params[b]
    vars = mapreduce(*, incident(bn, b, :call), init = p) do i
      j = bn[i, :arg]
      return symbolic_vars[j]
    end
  end

  infs = map(parts(bn, :Qout)) do tv
    flux = mapreduce(+, incident(bn, tv, :infusion), init = 0) do wa
      j = bn[wa, :influx]
      return ϕs[j]
    end
    flux -= mapreduce(+, incident(bn, tv, :effusion), init = 0) do wa
      j = bn[wa, :efflux]
      return ϕs[j]
    end
  end

  # We assume bnsir[:tanvar] ⊆ bnsir[:variable] here
  tanvar_idxs = indexin(bnsir[:tanvar], bnsir[:variable])
  zparts = zip(tanvar_idxs, infs)

  eqs = Equation[D(symbolic_vars[j::Int]) ~ rhs for (j, rhs) in zparts]
  ODESystem(eqs, t, symbolic_vars, symbolic_params, name=name)
end

# Compile bilayer network to a modeling toolkit expression
function compile(bn::Union{AbstractLabelledBilayerNetwork,AbstractBilayerNetwork})
  varstmt = :(@variables t)
  # get state names
  varnames = bnsir[:variable]
  # convert variables to time dependent variables and add to variable statements
  append!(varstmt.args, make_depvar.(bnsir[:variable], :t))

  # get transition names and add as parameters
  paramstmt = :(@parameters)
  params = bnsir[:parameter]
  append!(paramstmt.args, bnsir[:parameter])

  diffstmt = :(D = Differential(t))

  ϕs = map(parts(bn, :Box)) do b
    vars = map(incident(bn, b, :call)) do i
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
      return :(-$(Symbol("ϕ$j")))
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

#= /Users/fairbanksj/github/AlgebraicJulia/ASKEM-demos/MTK/petri.jl:39 =#
@variables t S(t) I(t) R(t)
@parameters inf rec
D = Differential(t)
ϕ1 = inf * S * I
ϕ2 = rec * I
eqs = [D(S) ~ +(-ϕ1), D(I) ~ ϕ1 + ϕ1 + -ϕ1 + -ϕ2, D(R) ~ +ϕ2]
example = ODESystem(eqs, t, name=:PetriNet)
@test ODESystem(bnsir) == example

# ... what is this
# (params...)
#
# begin
#   @variables t S(t) I(t) R(t)
#   @parameters inf rec
#   D = Differential(t)
#   ϕ1 = inf * S * I
#   ϕ2 = rec * I
#   eqs = [D(S) ~ +(-ϕ1), D(I) ~ ϕ1 + ϕ1 + -ϕ1 + -ϕ2, D(R) ~ +ϕ2]
#   ((S, I, R), (inf, rec), ODESystem(eqs, t, name=:PetriNet))
# end
