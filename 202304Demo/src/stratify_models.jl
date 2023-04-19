using AlgebraicPetri, AlgebraicPetri.TypedPetri
using Catlab.Programs, Catlab.Graphics
using Catlab.CategoricalAlgebra
using Catlab.WiringDiagrams

include("component_models.jl")

m1_uwd = @relation () where{(S::Pop, I::Pop, D::Pop, A::Pop, R::Pop, T::Pop, H::Pop, E::Pop)} begin
    infect(S,D,I,D)
    infect(S,A,I,A)    
    infect(S,R,I,R)
    infect(S,I,I,I)
    disease(I,D)
    disease(I,A)
    disease(I,H)
    disease(D,R)
    disease(D,H)
    disease(A,R)
    disease(A,H)
    disease(A,T)
    disease(R,T)
    disease(R,H)
    disease(T,E)
    disease(T,H)
end

m2_uwd = @relation () where{(S::Pop, E::Pop, I::Pop, A::Pop, H::Pop, R::Pop, D::Pop)} begin
    infect(S,I,E,I)
    infect(S,A,E,A)    
    infect(S,H,E,H)
    disease(E,I)
    disease(E,A)
    disease(I,H)
    disease(I,R)
    disease(I,D)
    disease(A,R)
    disease(A,D)
    disease(H,D)
    disease(H,R)
end

m3_uwd = @relation () where{(S::Pop, E::Pop, Iu::Pop, Ir::Pop, Q::Pop, R::Pop, D::Pop)} begin
    infect(S,Ir,E,Ir)
    infect(S,Iu,E,Iu)    
    infect(S,Ir,Q,Ir)
    infect(S,Iu,Q,Iu)
    disease(Q,Ir)
    disease(E,Ir)
    disease(E,Iu)
    disease(Ir,R)
    disease(Iu,R)
    disease(Q,S)
    disease(Ir,D)
end

oapply_mira_model(uwd) = oapply_typed(infectious_ontology, uwd, [Symbol("t$(n)") for n in 1:nboxes(uwd)])

m1_model = oapply_mira_model(m1_uwd)
m1_model = add_reflexives(
    m1_model,
    [[:strata], [:strata], [:strata], [:strata], [:strata], [], [:strata], []],
    infectious_ontology
)

m2_model = oapply_mira_model(m2_uwd)
m2_model = add_reflexives(
    m2_model,
    [[:strata], [:strata], [:strata], [:strata], [], [:strata], []],
    infectious_ontology
)

m3_model = oapply_mira_model(m3_uwd)
m3_model = add_reflexives(
    m3_model,
    [[:strata], [:strata], [:strata], [:strata], [], [:strata], []],
    infectious_ontology
)

disease_models = [m1_model, m2_model, m3_model]

policy_models = [nothing, vax_model, mask_model, mask_vax_model]

num_rgns = 2
base_travel_model = make_rgn(num_rgns)
travel_models = [nothing, base_travel_model, typed_product(base_travel_model, make_living(num_rgns))]

mdls = [typed_product(collect(filter(x -> !isnothing(x), pieces))) for pieces in Iterators.product(disease_models, policy_models, travel_models)]