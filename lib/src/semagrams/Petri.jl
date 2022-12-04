module Petri
export StratPetri, StratPetriSemagram, typed_petri

using ...Ontologies

using AlgebraicPetri
using Semagrams
using Catlab.Theories, Catlab.Present, Catlab.CSetDataStructures, Catlab.CategoricalAlgebra.CSets

"""
A petri net with extra features for stratification. Essentially, a Petri net typed by
the basic epidemiology ontology, but typed via attributes instead of via a morphism.

For writing down a Petri net for stratification in Semagrams.
"""
@present SchStratPetri <: SchLabelledReactionNet begin
  TransitionTypeValue::AttrType
  transitiontype::Attr(T,TransitionTypeValue)

  Stratification::AttrType
  stratificationwith::Attr(S,Stratification)
end

@acset_type StratPetriUntyped(SchStratPetri)

const StratPetri = StratPetriUntyped{Float64, Float64, Symbol, Symbol, Vector{Symbol}}

const petri_build = "https://semagrams-builds.s3.amazonaws.com/0d10593/petri/main.js"
const dev_petri_build = "http://localhost:8080/out/apps/petri/fullLinkJS.dest/main.js"

StratPetriSemagram(;dev=false) = Semagram{StratPetri}(
  dev ? dev_petri_build : petri_build,
  "Petri",
  Dict{Symbol,Function}(
    :Name => s -> Symbol(s),
    :Rate => s -> s,
    :Concentration => s -> s,
    :TransitionTypeValue => s ->
      if length(s) == 0
        :none
      else
        Symbol(s[1])
      end,
    :Stratification => s -> Symbol[Symbol(str) for str in s]
  ),
)

function get_transition_type(acs, i, transition_types)
  findfirst(transition_types .== subpart(acs, i, :transitiontype))
end

function typed_petri(sp::StratPetri, ontology)
  p = LabelledPetriNet()
  copy_parts!(p, sp)
  transition_types = Int[get_transition_type(sp, i, ontology[:tname]) for i in parts(sp, :T)]
  homomorphism(p, strip_names(ontology);
               initial=(T=transition_types,), type_components=(Name=x->nothing,))
end

end
