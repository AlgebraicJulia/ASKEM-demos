module EpidemiologicalPetriNets
export SIRD, SIRD_typed

using ..Ontologies

using AlgebraicPetri
using Catlab.CategoricalAlgebra

SIRD = LabelledPetriNet([:S, :I, :R, :D],
    :inf => ((:S,:I) => (:I,:I)),
    :recover => (:I=>:R),
    :death => (:I=>:D));

SIRD_typed = homomorphism(SIRD, strip_names(infectious_ontology);
                          initial=(T=[1,2,2],I=[1,2,3,3],O=[1,2,3,3]),
                          type_components=(Name=x->nothing,))

end
