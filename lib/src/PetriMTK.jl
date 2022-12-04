module PetriMTK
export rate_equation

using Catlab.CSetDataStructures
using AlgebraicPetri
using ModelingToolkit
using Symbolics: map_subscripts

function rate_equation(p::AbstractPetriNet; name)
  @parameters t
  sname′(i) = if has_subpart(p, :sname)
    sname(p, i)
  else
    Symbol("S", map_subscripts(i))
  end
  tname′(i) = if has_subpart(p, :tname)
    tname(p, i)
  else
    Symbol("r", map_subscripts(i))
  end
  S = [first(@variables $Si(t)) for Si in sname′.(1:ns(p))]
  r = [first(@variables $ri(t)) for ri in tname′.(1:nt(p))]
  D = Differential(t)

  tm = TransitionMatrices(p)

  coefficients = tm.output - tm.input

  transition_rates = [r[tr] * prod(S[s]^tm.input[tr,s] for s in 1:ns(p)) for tr in 1:nt(p)]

  concentration_eqs = [D(S[s]) ~ transition_rates' * coefficients[:,s] for s in 1:ns(p)]
  rate_eqs = [D(r[tr]) ~ 0.0 for tr in 1:nt(p)]
  eqs = [concentration_eqs; rate_eqs]

  ODESystem(eqs, t, [S; r], []; name)
end

# Next steps:
#
# The rates should be state variables too. Policies are then implemented via the event system.
# We can evaluate whether certain observables stay below a certain threshold using the event system
# as well.

end
