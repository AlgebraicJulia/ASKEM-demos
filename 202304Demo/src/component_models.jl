# ORIGINAL CODE WITH OLD API
# ===========================
# # Function to assemble a multi-region model
# function make_rgn(n,f_aug=false)
#     lstates = []
#     ltrans = []
#     for ii in 1:n
#         push!(lstates,Symbol("Rgn"*string(ii)))
#     end
#     for ii in 1:n
#         for jj in 1:n
#             if ii != jj
#                 push!(ltrans,Symbol("travel"*string(ii)*string(jj)) => ((lstates[ii])=>(lstates[jj])))
#             end
#         end
#     end
#     if f_aug
#         for ii in 1:n
#             push!(ltrans,Symbol("id_d"*string(ii)) => ((lstates[ii])=>(lstates[ii])))
#         end
#     end
#     if f_aug
#         for ii in 1:n
#             push!(ltrans,Symbol("id_i"*string(ii)) => ((lstates[ii],lstates[ii])=>(lstates[ii],lstates[ii])))
#         end
#     end
#     map(println,ltrans)
#     MultiRgn = LabelledPetriNet(lstates,ltrans...)  
# end
#
# # Multi-region augmented model
# MultiRgn_aug = make_rgn(num_rgns, true)
# AlgebraicPetri.Graph(MultiRgn_aug)
#
# # Typed MR model
# mr_S = repeat([s],ns(MultiRgn_aug))
# mr_T = repeat([t_strata],num_rgns*(num_rgns-1))
# mr_I = repeat([i_strata],num_rgns*(num_rgns-1))
# mr_O = repeat([o_strata],num_rgns*(num_rgns-1))
# for ii in 1:num_rgns
#     push!(mr_T,t_disease)
#     push!(mr_I,i_disease)
#     push!(mr_O,o_disease)
# end
# for ii in 1:num_rgns
#     push!(mr_T,t_interact)
#     push!(mr_I,i_interact1)
#     push!(mr_I,i_interact2)
#     push!(mr_O,o_interact1)
#     push!(mr_O,o_interact2)
# end
# MR_aug_typed = ACSetTransformation(MultiRgn_aug, types,
# 	  S = mr_S,
# 	  T = mr_T,
# 	  I = mr_I,
# 	  O = mr_O,
# 	  Name = name -> nothing # specify the mapping for the loose ACSet transform
# )
# @assert is_natural(MR_aug_typed)
# AlgebraicPetri.Graph(dom(MR_aug_typed))
#
# # Function to form stratified model
# typed_stratify(typed_model1, typed_model2) =
# 		Theories.compose(proj1(CategoricalAlgebra.pullback(typed_model1, typed_model2)), typed_model1)
#
# # Multi-region SIRD stratified model
# SIRD_MR = typed_stratify(SIRD_aug_typed, MR_aug_typed)
# AlgebraicPetri.Graph(dom(SIRD_MR))
#
#
# function make_living(n,f_aug=false)
#     lstates = []
#     ltrans = []
#     for ii in 1:n
#         push!(lstates,Symbol("Live"*string(ii)))
#     end
#     if f_aug
#         for ii in 1:n
#             push!(ltrans,Symbol("l_id_d"*string(ii)) => ((lstates[ii])=>(lstates[ii])))
#         end
#     end
#     if f_aug
#         for ii in 1:n
#             push!(ltrans,Symbol("l_id_s"*string(ii)) => ((lstates[ii])=>(lstates[ii])))
#         end
#     end
#     for ii in 1:n
#         for jj in 1:n
#             push!(ltrans,Symbol("l_interact"*string(ii)*string(jj)) => ((lstates[ii],lstates[jj])=>(lstates[ii],lstates[jj])))
#         end
#     end
#     map(println,ltrans)
#     MultiRgn = LabelledPetriNet(lstates,ltrans...)  
# end
#
# # Living augmented model
# Living_aug = make_living(num_rgns, true)
# AlgebraicPetri.Graph(Living_aug)
#
# # Typed Living model
# l_S = repeat([s],ns(Living_aug))
# l_T = repeat([t_disease],num_rgns)
# l_I = repeat([i_disease],num_rgns)
# l_O = repeat([o_disease],num_rgns)
# for ii in 1:num_rgns
#     push!(l_T,t_strata)
#     push!(l_I,i_strata)
#     push!(l_O,o_strata)
# end
# for ii in 1:num_rgns
#     for jj in 1:num_rgns
#         push!(l_T,t_interact)
#         push!(l_I,i_interact1)
#         push!(l_I,i_interact2)
#         push!(l_O,o_interact1)
#         push!(l_O,o_interact2)
#     end
# end
# Living_aug_typed = ACSetTransformation(Living_aug, types,
# 	  S = l_S,
# 	  T = l_T,
# 	  I = l_I,
# 	  O = l_O,
# 	  Name = name -> nothing # specify the mapping for the loose ACSet transform
# )
# @assert is_natural(Living_aug_typed)
# AlgebraicPetri.Graph(dom(Living_aug_typed))
#
# simple_trip_typed = typed_stratify(MR_aug_typed, Living_aug_typed)
# # AlgebraicPetri.Graph(dom(simple_trip_typed))
#    
# SIRD_trip = typed_stratify(SIRD_aug_typed,simple_trip_typed) 
# # AlgebraicPetri.Graph(dom(SIRD_trip))
#
#
#
# function make_nbr_rgn(ladj,f_aug=false)
#     n = length(ladj)
#     lstates = []
#     ltrans = []
#     for ii in 1:n
#         push!(lstates,Symbol("Rgn"*string(ii)))
#     end
#     if f_aug
#         for ii in 1:n
#             push!(ltrans,Symbol("id_s"*string(ii)) => ((lstates[ii])=>(lstates[ii])))
#         end
#     end
#     if f_aug
#         for ii in 1:n
#             push!(ltrans,Symbol("id_d"*string(ii)) => ((lstates[ii])=>(lstates[ii])))
#         end
#     end
#     for ii in 1:n
#         push!(ltrans,Symbol("i_local"*string(ii)) => ((lstates[ii],lstates[ii])=>(lstates[ii],lstates[ii])))
#         for jj in ladj[ii]
#             push!(ltrans,Symbol("i_travel"*string(ii)*string(jj)) => ((lstates[ii],lstates[jj])=>(lstates[ii],lstates[jj])))
#         end
#     end
#     map(println,ltrans)
#     nbr_rgn = LabelledPetriNet(lstates,ltrans...)  
# end
#
# # Multi-region augmented model
# num_rgns = 4
# NbrRgn_aug = make_nbr_rgn([[2,3],[1,4],[1,4],[2,3]], true)
# AlgebraicPetri.Graph(NbrRgn_aug)
#
# # Typed MR model
# mr_S = repeat([s],num_rgns)
# mr_T = repeat([t_strata],num_rgns)
# mr_I = repeat([i_strata],num_rgns)
# mr_O = repeat([o_strata],num_rgns)
# for ii in 1:num_rgns
#     push!(mr_T,t_disease)
#     push!(mr_I,i_disease)
#     push!(mr_O,o_disease)
# end
# for ii in 1:num_rgns
#     for jj in 1:(length(ladj[ii])+1)
#         push!(mr_T,t_interact)
#         push!(mr_I,i_interact1)
#         push!(mr_I,i_interact2)
#         push!(mr_O,o_interact1)
#         push!(mr_O,o_interact2)
#     end
# end
# Nbr_aug_typed = ACSetTransformation(NbrRgn_aug, types,
# 	  S = mr_S,
# 	  T = mr_T,
# 	  I = mr_I,
# 	  O = mr_O,
# 	  Name = name -> nothing # specify the mapping for the loose ACSet transform
# )
# @assert is_natural(Nbr_aug_typed)
# AlgebraicPetri.Graph(dom(Nbr_aug_typed))
#
# # Function to form stratified model
# typed_stratify(typed_model1, typed_model2) =
# 		Theories.compose(proj1(CategoricalAlgebra.pullback(typed_model1, typed_model2)), typed_model1)
#
# # Multi-region SIRD stratified model
# SIRD_MR = typed_stratify(SIRD_aug_typed, MR_aug_typed)
# AlgebraicPetri.Graph(dom(SIRD_MR))


# *******

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

Graph(infectious_ontology)

# Define Model

sir_uwd = @relation () where (S::Pop, I::Pop, R::Pop) begin
    infect(S, I, I, I)
    disease(I, R)
end

to_graphviz(sir_uwd, box_labels = :name, junction_labels = :variable)

tnames = [:beta, :gamma]
typed_sir = oapply_typed(infectious_ontology, sir_uwd, tnames)
Graph(dom(typed_sir))

typed_sir_aug = add_reflexives(
    typed_sir,
    [[:strata], [:strata], [:strata]],
    infectious_ontology
)

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

travel_2 = make_rgn(2)

Graph(dom(travel_2))

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

# Living augmented model

living_2 = make_living(2)

Graph(dom(living_2))

# Stratify SIR with travel model

num_rgns = 2
short_trip_model = typed_product([typed_sir_aug, make_rgn(num_rgns), make_living(num_rgns)])

Graph(dom(short_trip_model))
