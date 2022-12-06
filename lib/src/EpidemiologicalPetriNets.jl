module EpidemiologicalPetriNets
export SIR, SIRD, SEIR, SEIRD, SVIIvR,
  SIRD_typed

using ..Ontologies

using AlgebraicPetri
using Catlab.CategoricalAlgebra

SIR = LabelledPetriNet(
  [:S, :I, :R],
  :inf => ((:S,:I) => (:I,:I)),
  :recover => (:I=>:R)
)

SIRD = LabelledPetriNet(
  [:S, :I, :R, :D],
  :inf => ((:S,:I) => (:I,:I)),
  :recover => (:I=>:R),
  :death => (:I=>:D)
)

SEIR = LabelledPetriNet(
  [:S, :E, :I, :R],
  :inf => ((:S,:I) => (:E,:I)),
  :sicken => (:E=>:I),
  :recover => (:I=>:R),
)

SEIRD = LabelledPetriNet(
  [:S, :E, :I, :R, :D],
  :inf => ((:S,:I) => (:E,:I)),
  :sicken => (:E=>:I),
  :recover => (:I=>:R),
  :death => (:I=>:D)
)

SVIIvR = LabelledPetriNet(
  [:S, :V, :I, :Iv, :R],
  :vaccinate => (:S => :V),
  :si_inf => ((:S,:I) => (:I,:I)),
  :vi_inf => ((:V,:I) => (:Iv,:I)),
  :viv_inf => ((:V,:Iv) => (:Iv,:Iv)),
  :siv_inf => ((:S,:Iv) => (:I,:Iv)),
  :ivrecover => (:Iv=>:R),
  :irecover => (:I=>:R)
)

SIRD_typed = homomorphism(
  SIRD, strip_names(infectious_ontology);
  initial=(T=[1,2,2],I=[1,2,3,3],O=[1,2,3,3]),
  type_components=(Name=x->nothing,)
)

end
