module Stratify 
export stratify, stratify_except, age_strata, StrataSpec, add_cross_terms_with_rates

using Catlab.CategoricalAlgebra
using AlgebraicPetri

"""
Modify a typed petri net to add cross terms

NOTE this is just a helper function for `stratify`, below. It does not 
ever need to be called directly.
"""
function add_cross_terms(pn_crossterms, type_system)

  typed_pn, crossterms = deepcopy.(pn_crossterms)
  pn = dom(typed_pn)
  type_comps = Dict([k=>collect(v) for (k,v) in pairs(components(typed_pn))])
  for (s_i,cts) in enumerate(crossterms)
    for ct in cts 
      type_ind = findfirst(==(ct), type_system[:tname])
      is, os = [incident(type_system, type_ind, f) for f in [:it, :ot]]
      new_t = add_part!(pn, :T; tname=ct)
      add_parts!(pn, :I, length(is); is=s_i, it=new_t)
      add_parts!(pn, :O, length(os); os=s_i, ot=new_t)
      push!(type_comps[:T], type_ind)
      append!(type_comps[:I], is); append!(type_comps[:O], os); 
    end
  end
  return homomorphism(pn, codom(typed_pn); initial=type_comps, 
                      type_components=(Name=x->nothing,),)
end

"""
Perform a subset of possible petri net stratifications by indicating, for each 
species, which transitions in the *other* model you want that species to have.

For example, in stratifying a SIRD model with a Quarantine model, the args would 
be:
SIRD_typed=>[[:strata],[:strata],[:strata],[]],     # note: order is S I R D
Quarantine_typed=>[[:disease], [:disease,:infect]], # note: order is Q NQ 

This says that S,I,R can partake in strata-changing transitions (as defined in 
the Quarantine model) and Q,NQ can change disease state while only NQ can 
participate in infections. See the tests for more details.

This works by add cross terms prior to taking a pullback.
"""
function stratify(pn1, pn2, type_system)
  pullback([add_cross_terms(pn, type_system) for pn in [pn1, pn2]]) |> apex
end

function stratify_except(pn1, pn2, type_system)
  cross1, cross2 = map([pn1=>pn2,pn2=>pn1]) do (pn, other) 
    map(parts(dom(pn[1]),:S)) do s 
      setdiff(type_system[other[1][:T]|>collect,:tname], get(Dict(pn[2]), s, []))
    end
  end
  stratify(pn1[1]=>cross1,pn2[1]=>cross2, type_system)
end


#=
THESE ARE STRATIFY FUNCTIONS FROM Structured-Epidemic-Modeling
stratify(typed_model1, typed_model2) = ob(pullback(typed_model1, typed_model2))

typed_stratify(typed_model1, typed_model2) =
  compose(proj1(pullback(typed_model1, typed_model2)), typed_model1);
=#

abstract type AbstractStrataSpec end

struct StrataSpec <: AbstractStrataSpec
  tpn::ACSetTransformation
  tlist::Vector{Vector{Symbol}}
end

function add_cross_terms_with_rates(ss::StrataSpec, type_system)
  typed_pn = deepcopy(ss.tpn)
  crossterms = deepcopy(ss.tlist)
  pn = dom(typed_pn)
  type_comps = Dict([k=>collect(v) for (k,v) in pairs(components(typed_pn))])
  for (s_i,cts) in enumerate(crossterms)
    for ct in cts 
      type_ind = findfirst(==(ct), type_system[:tname])
      is, os = [incident(type_system, type_ind, f) for f in [:it, :ot]]
      new_t = isa(pn, LabelledReactionNet) ? add_part!(pn, :T; tname=ct, rate=1.0) : add_part!(pn, :T; tname=ct)
      add_parts!(pn, :I, length(is); is=s_i, it=new_t)
      add_parts!(pn, :O, length(os); os=s_i, ot=new_t)
      push!(type_comps[:T], type_ind)
      append!(type_comps[:I], is); append!(type_comps[:O], os); 
    end
  end
  return homomorphism(pn, codom(typed_pn); initial=type_comps, 
                        type_components=(Name=x->nothing,Rate=x->nothing, Concentration=x->nothing),)
end


"""
Create a petri net with `n` species, a unary transition on each, and a binary
transition for each (unordered) pair of species.
"""
function age_strata(n::Int)
  res = oplus(fill(terminal(PetriNet) |> apex, n))
  for ij in filter(ij-> ij[1] <= ij[2], collect(Iterators.product(1:n, 1:n)))
    t = add_part!(res, :T)
    add_parts!(res, :I, 2; it=t, is=ij); add_parts!(res, :O, 2; ot=t, os=ij)
  end
  return res
end 


end # module 
