#=module Oct2022Demo

export stratify, typed_stratify, simulate, calibrate, MakeReactionSystem, generate_data, get_infected_states, get_susceptible_states,
    get_recovered_states, optimise_p, full_train, obs_SIRD, makeLossFunction, plot_obs_w_ests, full_train_rand
 =#

using AlgebraicPetri
using Catlab
using Catlab.CategoricalAlgebra
using JSON
using CSV

#=
THESE ARE STRATIFY FUNCTIONS FROM Structured-Epidemic-Modeling
stratify(typed_model1, typed_model2) = ob(pullback(typed_model1, typed_model2))

typed_stratify(typed_model1, typed_model2) =
  compose(proj1(pullback(typed_model1, typed_model2)), typed_model1);
=#

abstract type AbstractStrataSpec end

struct StrataSpec <: AbstractStrataSpec
  tpn::ACSetTransformation
  tlist::Vector{Vector{Symbol}}
end

"""Modify a typed petri net to add cross terms"""
#=function add_cross_terms(pn_crossterms, type_system)

  typed_pn, crossterms = deepcopy.(pn_crossterms)
  pn = dom(typed_pn)
  type_comps = Dict([k=>collect(v) for (k,v) in pairs(components(typed_pn))])
  for (s_i,cts) in enumerate(crossterms)
    for ct in cts 
      type_ind = findfirst(==(ct), type_system[:tname])
      is, os = [incident(type_system, type_ind, f) for f in [:it, :ot]]
      new_t = add_part!(pn, :T; tname=ct)
      add_parts!(pn, :I, length(is); is=s_i, it=new_t)
      add_parts!(pn, :O, length(os); os=s_i, ot=new_t)
      push!(type_comps[:T], type_ind)
      append!(type_comps[:I], is); append!(type_comps[:O], os); 
    end
  end
  return homomorphism(pn, codom(typed_pn); initial=type_comps, 
                        type_components=(Name=x->nothing,),)
end

function add_cross_terms(ss::StrataSpec, type_system)
  return add_cross_terms(ss.tpn=>ss.tlist, type_system)
end=#

function add_cross_terms!(ss::StrataSpec, type_system)
  typed_pn = ss.tpn
  crossterms = ss.tlist
  pn = dom(typed_pn)
  type_comps = Dict([k=>collect(v) for (k,v) in pairs(components(typed_pn))])
  for (s_i,cts) in enumerate(crossterms)
    for ct in cts 
      type_ind = findfirst(==(ct), type_system[:tname])
      is, os = [incident(type_system, type_ind, f) for f in [:it, :ot]]
      new_t = add_part!(pn, :T; tname=ct)
      add_parts!(pn, :I, length(is); is=s_i, it=new_t)
      add_parts!(pn, :O, length(os); os=s_i, ot=new_t)
      push!(type_comps[:T], type_ind)
      append!(type_comps[:I], is); append!(type_comps[:O], os); 
    end
  end
  return homomorphism(pn, codom(typed_pn); initial=type_comps, 
                        type_components=(Name=x->nothing,),)
end

function add_cross_terms(ss::StrataSpec, type_system)
  typed_pn = deepcopy(ss.tpn)
  crossterms = deepcopy(ss.tlist)
  pn = dom(typed_pn)
  type_comps = Dict([k=>collect(v) for (k,v) in pairs(components(typed_pn))])
  for (s_i,cts) in enumerate(crossterms)
    for ct in cts 
      type_ind = findfirst(==(ct), type_system[:tname])
      is, os = [incident(type_system, type_ind, f) for f in [:it, :ot]]
      new_t = add_part!(pn, :T; tname=ct)
      add_parts!(pn, :I, length(is); is=s_i, it=new_t)
      add_parts!(pn, :O, length(os); os=s_i, ot=new_t)
      push!(type_comps[:T], type_ind)
      append!(type_comps[:I], is); append!(type_comps[:O], os); 
    end
  end
  return homomorphism(pn, codom(typed_pn); initial=type_comps, 
                        type_components=(Name=x->nothing,),)
end
  
function add_cross_terms_with_rates(ss::StrataSpec, type_system)
  typed_pn = deepcopy(ss.tpn)
  crossterms = deepcopy(ss.tlist)
  pn = dom(typed_pn)
  type_comps = Dict([k=>collect(v) for (k,v) in pairs(components(typed_pn))])
  for (s_i,cts) in enumerate(crossterms)
    for ct in cts 
      type_ind = findfirst(==(ct), type_system[:tname])
      is, os = [incident(type_system, type_ind, f) for f in [:it, :ot]]
      new_t = isa(pn, LabelledReactionNet) ? add_part!(pn, :T; tname=ct, rate=1.0) : add_part!(pn, :T; tname=ct)
      add_parts!(pn, :I, length(is); is=s_i, it=new_t)
      add_parts!(pn, :O, length(os); os=s_i, ot=new_t)
      push!(type_comps[:T], type_ind)
      append!(type_comps[:I], is); append!(type_comps[:O], os); 
    end
  end
  return homomorphism(pn, codom(typed_pn); initial=type_comps, 
                        type_components=(Name=x->nothing,Rate=x->nothing, Concentration=x->nothing),)
end

"""Add cross terms before taking pullback"""
function stratify_span(ss1, ss2, type_system)
  #type_system = codom(ss1.tpn)
  pb = Catlab.CategoricalAlgebra.pullback(add_cross_terms(ss1,type_system), add_cross_terms(ss2, type_system))
  return pb
end

stratify(ss1, ss2, type_system) = begin
  pb = stratify_span(ss1,ss2,type_system)
  f = proj1(pb)
  return apex(pb), makeObsFunction(f)
end

counter(a) = [count(==(i),a) for i in unique(a)]

@present SchStratPetri <: TheoryLabelledReactionNet begin
  TransitionTypeValue::AttrType
  transitiontype::Attr(T,TransitionTypeValue)

  Stratification::AttrType
  stratificationwith::Attr(S,Stratification)
end

@acset_type StratPetriUntyped(SchStratPetri)

const StratPetri = StratPetriUntyped{Float64, Float64, Symbol, Symbol, Vector{Symbol}}

strat_petri_decoders = Dict{Symbol,Function}(
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
)

function gettype(s::StratPetri, i)
  res = findfirst([:strata, :disease, :type] .== subpart(s, i, :transitiontype))
  if !isnothing(res)
    res
  else
    0
  end
end

function tostrataspec(s::StratPetri, t)
  p = LabelledPetriNet()
  copy_parts!(p, s)
  types = Int[gettype(s, i) for i in parts(s, :T)]
  p_typed = homomorphism(
    p, t;
    initial=(T=types,), type_components=(Name=x->nothing,)
  )
  strat_p = StrataSpec(p_typed, subpart(s, :stratificationwith))
end

function load_mira(fname)
  mdl_orig_mira = read_json_acset(OrigMIRANet{Any,Any,Any,Any,Any,Any}, fname)
  mdl_orig_mira[:sname] .= Symbol.(mdl_orig_mira[:sname])
  mdl_orig_mira[:tname] .= Symbol.(mdl_orig_mira[:tname])
  mdl_orig = LabelledPetriNet(mdl_orig_mira);
  return mdl_orig
end

function load_mira_curated(fname, fillrate=1.0)
  mdl_disease_mira = read_json_acset(MIRANet{Any,Any,Any,Any,Any,Any}, fname)
  input_rates = deepcopy(mdl_disease_mira[:rate])
	min_rate = minimum(input_rates[map(!isnothing,input_rates)])
  filled_rates = input_rates
  filled_rates[map(isnothing,input_rates)] = repeat([0.1],sum(map(isnothing,input_rates)))
  mdl_disease_mira[:rate] = filled_rates
  mdl_disease_mira[:concentration] = 0.0
  mdl_disease_mira[:sname] .= Symbol.(mdl_disease_mira[:sname])
  mdl_disease_mira[:tname] .= Symbol.(mdl_disease_mira[:tname])
  return mdl_disease_mira
end

get_infected_states(g::AbstractLabelledPetriNet) =
   [i for (i,s) in enumerate(g[:sname]) if occursin("I", string(s))]

get_susceptible_states(g::AbstractLabelledPetriNet) =
   [i for (i,s) in enumerate(g[:sname]) if occursin("S", string(s))]

get_recovered_states(g::AbstractLabelledPetriNet) =
   [i for (i,s) in enumerate(g[:sname]) if occursin("R", string(s)) || occursin("V", string(s))]

get_dead_states(g::AbstractLabelledPetriNet) =
   [i for (i,s) in enumerate(g[:sname]) if occursin("D", string(s))]

#=function obsSIRD(model::AbstractLabelledPetriNet, sol, sample_times)
    inf_sample_vals = sol(sample_times)[get_infected_states(model),:]
    susc_sample_vals = sol(sample_times)[get_susceptible_states(model),:]
    rec_sample_vals = sol(sample_times)[get_recovered_states(model),:]
    dead_sample_vals = sol(sample_times)[get_dead_states(model),:]

    total_inf_samples = map(sum, eachcol(inf_sample_vals))
    total_susc_samples = map(sum, eachcol(susc_sample_vals))
    total_rec_samples = map(sum, eachcol(rec_sample_vals))
    total_dead_samples = map(sum, eachcol(dead_sample_vals))

    labels = reshape(["I", "R or V", "S", "D"],1,4)

    return hcat(total_inf_samples, total_rec_samples, total_susc_samples, total_dead_samples), labels
end=#

function state_preimage(f::ACSetTransformation, y)
  fS = f.components.S
  return [ x for x in dom(fS) if fS(x) == y ]
end

function pushforwardObs(f::ACSetTransformation,sol,sample_times)

  state_samples = sol(sample_times)
  obs_samples = []
  for ii in 1:ns(codom(f))
    tmp = sol(sample_times)[state_preimage(f, ii),:]
    push!(obs_samples,map(sum, eachcol(tmp)))
  end

  return hcat(obs_samples...), snames(codom(f))
end

function makeObsFunction(f::ACSetTransformation)
  return (sol, sample_times) -> pushforwardObs(f, sol, sample_times)
end
