module Stratify 
export stratify, stratify_except,stratify_typed,stratify_except_typed, 
       age_strata, StrataSpec, add_cross_terms_with_rates, decompose

using Catlab.CategoricalAlgebra
using AlgebraicPetri
using ..Ontologies: strip_names

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

Returns a typed model, i.e. a map in Petri.
"""
function stratify_typed(pn1, pn2, type_system)
  pn1′, pn2′ = [add_cross_terms(pn, type_system) for pn in [pn1, pn2]]
  pb = pullback(pn1′, pn2′) 
  return first(legs(pb)) ⋅ pn1′
end

function stratify_except_typed(pn1, pn2, type_system)
  cross1, cross2 = map([pn1=>pn2,pn2=>pn1]) do (pn, other) 
    map(parts(dom(pn[1]),:S)) do s 
      setdiff(type_system[other[1][:T]|>collect,:tname], get(Dict(pn[2]), s, []))
    end
  end
  stratify_typed(pn1[1]=>cross1,pn2[1]=>cross2, type_system) 
end

stratify(pn1, pn2, type_system) = stratify_typed(pn1,pn2,type_system) |> dom

stratify_except(pn1, pn2, type_system) = 
  stratify_except_typed(pn1,pn2,type_system) |> dom



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



"""Copied from Catlab PR 710"""
function migrate(f::TightACSetTransformation, F::DeltaMigration) 
  d = Dict()
  for (ob_dom,ob_codom) in F.functor.ob_map
    if Symbol(ob_codom) ∈ keys(components(f))
      d[Symbol(ob_dom)] = f[Symbol(ob_codom)]
    end
  end
  TightACSetTransformation(NamedTuple(d), F(dom(f)), F(codom(f)))
end

"""Remove attributes from a LabelledPetriNet"""
strip_names_F = DeltaMigration(FinFunctor(
  Dict(:S=>:S,:T=>:T,:O=>:O,:I=>:I), 
  Dict(k => k for k in SchPetriNet.generators[:Hom]), 
  SchPetriNet, SchLabelledPetriNet), LabelledPetriNet, PetriNet)

"""Convert a LabeledPetriNet morphism into a PetriNet morphism, ignoring attrs"""
function remove_names(p::ACSetTransformation)
  init = NamedTuple([k=>collect(v) for (k,v) in pairs(components(p))])
  dom_codom = strip_names.([dom(p), codom(p)]) # turn into TIGHT transformations
  migrate(homomorphism(dom_codom...; initial=init), strip_names_F)
end


"""
Given a target model, determine if it is stratified by any pair of 
models from a list of candidates. 
"""
function decompose(tgt_model, candidate_legs) 
  stripped_tgt_model = remove_names(tgt_model)
  n = length(candidate_legs)
  res = []
  for i in 1:n 
    icand = remove_names(candidate_legs[i])
    for j in i:n 
      jcand = remove_names(candidate_legs[j])
      pb = first(legs(pullback(icand, jcand))) ⋅ icand
      if any(isomorphisms(dom(stripped_tgt_model), dom(pb))) do iso # copied from PR 709
          force(stripped_tgt_model) == force(iso ⋅ pb)
        end
          push!(res, candidate_legs[i]=>candidate_legs[j])
      end
    end
  end 
  return res
end

end # module 
