using ASKEM
using Catlab, Catlab.CategoricalAlgebra
using AlgebraicPetri
using Test


"""Copied from PR 710"""
function migrate(f::TightACSetTransformation, F::DeltaMigration) 
  d = Dict()
  for (ob_dom,ob_codom) in F.functor.ob_map
    if Symbol(ob_codom) ∈ keys(components(f))
      d[Symbol(ob_dom)] = f[Symbol(ob_codom)]
    end
  end
  TightACSetTransformation(NamedTuple(d), F(dom(f)), F(codom(f)))
end

strip_names_F = DeltaMigration(FinFunctor(
  Dict(:S=>:S,:T=>:T,:O=>:O,:I=>:I), 
  Dict(k => k for k in SchPetriNet.generators[:Hom]), 
  SchPetriNet, SchLabelledPetriNet), LabelledPetriNet, PetriNet)

"""Convert a LabeledPetriNet morphism into a PetriNet morphism, ignoring attrs"""
function remove_names(p::ACSetTransformation)
  init = NamedTuple([k=>collect(v) for (k,v) in pairs(components(p))])
  dom_codom = strip_names.([dom(p), codom(p)]) # turn into TIGHT transformations
  migrate(homomorphism(dom_codom...; initial=init), 
          strip_names_F)
end


"""
Given a target model, we wish to determine if it is stratified by any model of 
a list of candidates. We assume we have one of the two legs of the pullback. 
"""
function decompose(tgt_model, leg_1, candidate_leg2s) 
  stripped_tgt_model, stripped_leg_1 = remove_names.([tgt_model, leg_1])
  filter(candidate_leg2s) do cand
    pb = first(legs(pullback(stripped_leg_1, remove_names(cand)))) ⋅ stripped_leg_1
    any(isomorphisms(dom(stripped_tgt_model), dom(pb))) do iso # copied from PR 709
      force(stripped_tgt_model) == force(iso ⋅ pb)
    end
  end 
end


# Example
#########

# Two candidates for disease dynamics
SIR = LabelledPetriNet([:S, :I, :R],
    :inf => ((:S,:I) => (:I,:I)),
    :recover => (:I=>:R),
    :id => (:S=>:S), :id=>(:I=>:I),:id=>(:R=>:R))

SIRD = LabelledPetriNet([:S, :I, :R, :D],
    :inf => ((:S,:I) => (:I,:I)),
    :recover => (:I=>:R),
    :death => (:I=>:D),
    :id => (:S=>:S), :id=>(:I=>:I),:id=>(:R=>:R))

# One possible strata considered
Quarantine = LabelledPetriNet([:Q,:NQ],
    :inf => ((:NQ,:NQ)=>(:NQ,:NQ)),
    :id=>((:Q)=>(:Q)), 
    :id => ((:NQ) => (:NQ)),
    :quarantine => ((:NQ)=>(:Q)),
    :unquarantine => ((:Q)=>(:NQ)))

# Add type information
SIR_typed = homomorphism(SIR, strip_names(infectious_ontology);
    initial=(T=[1,2,3,3,3],I=[1,2,3,4,4,4],O=[1,2,3,4,4,4]),
    type_components=(Name=x->nothing,))

SIRD_typed = homomorphism(SIRD, strip_names(infectious_ontology);
    initial=(T=[1,2,2,3,3,3],I=[1,2,3,3,4,4,4],O=[1,2,3,3,4,4,4]),
    type_components=(Name=x->nothing,))

Quarantine_typed = homomorphism(Quarantine, strip_names(infectious_ontology);
    initial=(T=[1,2,2,3,3],I=Dict(1=>1,2=>2),O=Dict(1=>1,2=>2)), 
    type_components=(Name=x->nothing,))

# The "true" model is SIRD x Quarantine
tgt_model = first(legs(pullback(SIRD_typed, Quarantine_typed))) ⋅ SIRD_typed

# We give both disease dynamics as possibilities but only the correct one is returned
@test only(decompose(tgt_model, Quarantine_typed, [SIR_typed, SIRD_typed])) == SIRD_typed

# Stratification order doesn't matter 
tgt_model = first(legs(pullback(Quarantine_typed,SIRD_typed))) ⋅ Quarantine_typed
@test only(decompose(tgt_model, Quarantine_typed, [SIR_typed, SIRD_typed])) == SIRD_typed
