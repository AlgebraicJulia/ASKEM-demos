using ASKEM
using Catlab, Catlab.CategoricalAlgebra
using AlgebraicPetri
using Test

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

"""An isomorphic version of SIRD that has different names and order of parts."""
SIRD2 = LabelledPetriNet([:Sus, :Inf, :Dea, :Rec],
    :d => (:Inf=>:Dea),
    :i => ((:Sus,:Inf) => (:Inf,:Inf)),
    :r => (:Inf=>:Rec),
    :id => (:Sus=>:Sus), :id=>(:Inf=>:Inf),:id=>(:Rec=>:Rec))


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

SIRD2_typed = homomorphism(SIRD2, strip_names(infectious_ontology);
    initial=(T=[2,1,2,3,3,3],I=[3,1,2,3,4,4,4],O=[3,1,2,3,4,4,4]),
    type_components=(Name=x->nothing,))

Quarantine_typed = homomorphism(Quarantine, strip_names(infectious_ontology);
    initial=(T=[1,2,2,3,3],I=Dict(1=>1,2=>2),O=Dict(1=>1,2=>2)), 
    type_components=(Name=x->nothing,))

# The "true" model is SIRD x Quarantine
tgt_model = first(legs(pullback(SIRD_typed, Quarantine_typed))) ⋅ SIRD_typed

res = decompose(tgt_model, [SIR_typed, SIRD_typed, Quarantine_typed, SIRD2_typed])

# order doesn't matter, naming doesn't matter
@test res == [(SIRD_typed => Quarantine_typed),
              (Quarantine_typed => SIRD2_typed)]