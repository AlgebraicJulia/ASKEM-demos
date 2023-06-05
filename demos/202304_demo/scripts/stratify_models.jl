using DrWatson
@quickactivate "202304_demo"

println(
"""
Currently active project is: $(projectname())

Path of active project: $(projectdir())
"""
)

using AlgebraicPetri, AlgebraicPetri.TypedPetri
using Catlab.Programs, Catlab.Graphics
using Catlab.CategoricalAlgebra
using Catlab.WiringDiagrams

include(srcdir("component_models.jl"))

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

###
# Compute runtimes and number of states and transitions as a function of number of regions for models with simple trip
###

using BenchmarkTools
num_states = []
num_trans = []
num_states_3way = []
num_trans_3way = []
run_times = []
run_times_3way = []
l_num_rgns = [2, 4, 8, 16, 32, 50]
for num_rgns in l_num_rgns
    println(num_rgns)
    tmp_travel_model = make_rgn(num_rgns)
    tmp_simple_trip = typed_product(tmp_travel_model, make_living(num_rgns))
    println(ns(dom(tmp_simple_trip)))
    # push!(run_times, @elapsed typed_product(disease_models[1],tmp_simple_trip))
    tmp_time = @elapsed tmp_strat = typed_product(disease_models[1],tmp_simple_trip)
    push!(run_times, tmp_time)
    push!(num_states, ns(dom(tmp_strat)))
    push!(num_trans, nt(dom(tmp_strat)))
    # push!(run_times_3way, @elapsed typed_product([disease_models[1],policy_models[4],tmp_simple_trip]))
    tmp_time = @elapsed tmp_strat = typed_product([disease_models[1],policy_models[4],tmp_simple_trip])
    push!(run_times_3way, tmp_time)
    push!(num_states_3way, ns(dom(tmp_strat)))
    push!(num_trans_3way, nt(dom(tmp_strat)))
end

using LsqFit
@. quadmodel(x, p) = p[1]+p[2]*x + p[3]x^[2]
p0 = [0,.5,.5]
quadfit1_1 = LsqFit.curve_fit(quadmodel, l_num_rgns[2:end], run_times[2:end], p0)
quadfit1_2 = LsqFit.curve_fit(quadmodel, l_num_rgns, run_times_3way, p0)

@. quadmodel2(x, p) = p[1]+p[2]x^[2]
p0 = [0,.5]
quadfit2_1 = LsqFit.curve_fit(quadmodel2, l_num_rgns[2:end], run_times[2:end], p0)
quadfit2_2 = LsqFit.curve_fit(quadmodel2, l_num_rgns, run_times_3way, p0)

using Plots
plot(l_num_rgns,run_times;linewidth=2,marker=(:circle,4),labels="SIDARTHE + Simple Trip")
plot!(l_num_rgns,run_times_3way;linewidth=2,marker=(:circle,4),labels="SIDARTHE + MaskVax + Simple Trip")
plot!(title = "Runtime vs Number of Geographic Regions", xlabel = "Number of Geographic Regions", ylabel = "Runtime (sec)")
savefig(plotsdir("runtime_vs_num_rgns.svg"))

plot(l_num_rgns,num_states;linewidth=2,marker=(:circle,4),labels="SIDARTHE + Simple Trip")
plot!(l_num_rgns,num_states_3way;linewidth=2,marker=(:circle,4),labels="SIDARTHE + MaskVax + Simple Trip")
plot!(title = "Num of States vs Num of Geographic Regions", xlabel = "Number of Geographic Regions", ylabel = "Number of States")
savefig(plotsdir("num_states_vs_num_rgns.svg"))

plot(l_num_rgns,num_trans;linewidth=2,marker=(:circle,4),labels="SIDARTHE + Simple Trip")
plot!(l_num_rgns,num_trans_3way;linewidth=2,marker=(:circle,4),labels="SIDARTHE + MaskVax + Simple Trip")
plot!(title = "Num of Transitions vs Num of Geographic Regions", xlabel = "Number of Geographic Regions", ylabel = "Number of Transitions")
savefig(plotsdir("num_trans_vs_num_rgns.svg"))
