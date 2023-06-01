using AlgebraicPetri, AlgebraicPetri.TypedPetri
using Catlab.Programs, Catlab.Graphics
using Catlab.CategoricalAlgebra

const infectious_ontology = LabelledPetriNet(
    [:Pop],
    :infect => ((:Pop, :Pop) => (:Pop, :Pop)),
    :disease => (:Pop => :Pop),
    :strata => (:Pop => :Pop)
)

Graph(infectious_ontology)

sir_uwd = @relation () where (S::Pop, I::Pop, R::Pop) begin
    infect(S, I, I, I)
    disease(I, R)
end

to_graphviz(sir_uwd, box_labels = :name, junction_labels = :variable)

tnames = [:beta, :gamma]
typed_sir = oapply_typed(infectious_ontology, sir_uwd, tnames)
Graph(dom(typed_sir))

N = 2
snames = [Symbol("Age$i") for i in 1:N]

typed_age = pairwise_id_typed_petri(infectious_ontology, :Pop, :infect, snames)

Graph(dom(typed_age))

typed_age_aug = add_reflexives(
    typed_age,
    repeat([[:disease, :strata]], N),
    infectious_ontology
)

Graph(dom(typed_age_aug))

typed_sir_aug = add_reflexives(
    typed_sir,
    [[:strata], [:strata], [:strata]],
    infectious_ontology
)

Graph(dom(typed_sir_aug))

M = 2
snames = [Symbol("Country$i") for i in 1:M]

# Need a function like `pairwise_id_typed_petri` but without identity
travel_uwd = @relation () where (Country1::Pop, Country2::Pop) begin
    strata(Country1, Country2)
    strata(Country2, Country1)
end

to_graphviz(travel_uwd, box_labels = :name, junction_labels = :variable)

typed_travel = oapply_typed(infectious_ontology, travel_uwd, [:Country1_Country2, :Country2_Country1])

Graph(dom(typed_travel))

# Augment the travel model with reflexives

typed_travel_aug = add_reflexives(
    typed_travel,
    repeat([[:infect, :disease]], M),
    infectious_ontology
)

# Graph(dom(typed_sir_age_aug))

typed_sir_age_travel = typed_product([typed_sir_aug, typed_age_aug, typed_travel_aug])

Graph(dom(typed_sir_age_travel))
