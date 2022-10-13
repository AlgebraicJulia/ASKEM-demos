using Catlab.CategoricalAlgebra
using Catlab.Present, Catlab.Theories
using AlgebraicPetri

@present SchStratPetri <: TheoryLabelledReactionNet begin
  TransitionTypeValue::AttrType
  transitiontype::Attr(T,TransitionTypeValue)

  Stratification::AttrType
  stratificationwith::Attr(S,Stratification)
end

@acset_type StratPetriUntyped(SchStratPetri)

const StratPetri = StratPetriUntyped{Float64, Float64, Symbol, Symbol, Vector{Symbol}}
