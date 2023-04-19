using AlgebraicPetri, AlgebraicPetri.TypedPetri
using Catlab.Programs, Catlab.Graphics
using Catlab.CategoricalAlgebra
using Catlab.WiringDiagrams

# Define Type System

const infectious_ontology = LabelledPetriNet(
  [:Pop],
  :infect=>((:Pop, :Pop)=>(:Pop, :Pop)),
  :disease=>(:Pop=>:Pop),
  :strata=>(:Pop=>:Pop)
)

# # Define Model

# sir_uwd = @relation () where (S::Pop, I::Pop, R::Pop) begin
#     infect(S, I, I, I)
#     disease(I, R)
# end

# to_graphviz(sir_uwd, box_labels = :name, junction_labels = :variable)

# tnames = [:beta, :gamma]
# typed_sir = oapply_typed(infectious_ontology, sir_uwd, tnames)
# Graph(dom(typed_sir))

# typed_sir_aug = add_reflexives(
#     typed_sir,
#     [[:strata], [:strata], [:strata]],
#     infectious_ontology
# )

# Travel Model

function make_rgn(n)
    uwd = RelationalPrograms.TypedUnnamedRelationDiagram{Symbol,Symbol,Symbol}()
    # Add junctions and store their mapped keys
    junctions = Dict(begin
        junction = Symbol("Region$(i)")
        junction => add_junction!(uwd, :Pop, variable=junction)
    end for i in 1:n)

    # figure out unique pairs of junctions
    pairs = filter(x -> first(x) != last(x), collect(Iterators.product(keys(junctions), keys(junctions))))
    for pair in pairs
        box = add_box!(uwd, [junction_type(uwd, junctions[p]) for p in pair], name=:strata)
        for (rgn, port) in zip(pair, ports(uwd, box))
            set_junction!(uwd, port, junctions[rgn])
        end
    end
    # convert uwd to ACSetTransformation using type ontology
    act = oapply_typed(infectious_ontology, uwd, [Symbol("$(a)_$(b)") for (a,b) in pairs])
    # if reflexives are provided, add them `n` times
    add_reflexives(act, repeat([[:infect,:disease]], n), infectious_ontology)
end

# Living Model

function make_living(n)
    states = [Symbol("Living$(i)") for i in 1:n]
    typed_living = pairwise_id_typed_petri(infectious_ontology, :Pop, :infect, states)
    add_reflexives(
        typed_living,
        repeat([[:disease, :strata]], n),
        infectious_ontology
    )
end

# Masking model

masking_uwd = @relation () where (M::Pop, UM::Pop) begin
    disease(M, UM)    
    disease(UM, M)
    infect(M, UM, M, UM)    
    infect(UM, UM, UM, UM)
end
tnames = [:unmask, :mask, :infect_um, :infect_uu]
mask_model = oapply_typed(infectious_ontology, masking_uwd, tnames)

mask_model = add_reflexives(
    mask_model,
    [[:strata], [:strata]],
    infectious_ontology
)

# Vaccine model

vax_uwd = @relation () where (UV::Pop, V::Pop) begin
    strata(UV, V)
    infect(V, V, V, V)    
    infect(V, UV, V, UV)    
    infect(UV, V, UV, V)    
    infect(UV, UV, UV, UV)
end
tnames = [:vax, :infect_vv, :infect_uv, :infect_vu, :infect_uu]
vax_model = oapply_typed(infectious_ontology, vax_uwd, tnames)

vax_model = add_reflexives(
    vax_model,
    [[:disease], [:disease]],
    infectious_ontology
)


# Mask-Vax Model

mask_vax_uwd = @relation () where (UV_UM::Pop, UV_M::Pop, V_UM::Pop, V_M::Pop) begin
    strata(UV_UM, UV_M)
    strata(UV_M, UV_UM)
    strata(V_UM, V_M)
    strata(V_M, V_UM)
    strata(UV_UM,V_UM)
    strata(UV_M,V_M)
    infect(V_UM, V_UM, V_UM, V_UM)    
    infect(V_UM, UV_UM, V_UM, UV_UM)    
    infect(UV_UM, V_UM, UV_UM, V_UM)    
    infect(UV_UM, UV_UM, UV_UM, UV_UM)
    infect(V_M, V_UM, V_M, V_UM)    
    infect(V_M, UV_UM, V_M, UV_UM)    
    infect(UV_M, V_UM, UV_M, V_UM)    
    infect(UV_M, UV_UM, UV_M, UV_UM)
end
tnames = [:mask_uv, :unmask_uv, :mask_v, :unmask_v, :vax_um, :vax_m, :infect_um_vv, :infect_um_uv, :infect_um_vu, :infect_um_uu, :infect_m_vv, :infect_m_uv, :infect_m_vu, :infect_m_uu]
mask_vax_model = oapply_typed(infectious_ontology, mask_vax_uwd, tnames)

mask_vax_model = add_reflexives(
    mask_vax_model,
    [[:disease], [:disease], [:disease], [:disease]],
    infectious_ontology
)
