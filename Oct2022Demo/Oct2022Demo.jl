#=module Oct2022Demo

export stratify, typed_stratify, simulate, calibrate, MakeReactionSystem, generate_data, get_infected_states, get_susceptible_states,
    get_recovered_states, optimise_p, full_train, obs_SIRD, makeLossFunction, plot_obs_w_ests, full_train_rand
 =#

using AlgebraicPetri
using Catlab
using Catlab.CategoricalAlgebra
using Catalyst
using OrdinaryDiffEq
using DiffEqFlux
using Plots
using JSON
using CSV
using DataFrames

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
function MakeReactionSystem(pn::AbstractPetriNet)
    @parameters t k[1:(nt(pn) + ns(pn))]
    @variables S[1:ns(pn)](t)

    rxs = map(1:nt(pn)) do t
      inpts = pn[incident(pn, t, :it),:is]
      otpts = pn[incident(pn, t, :ot),:os]
      in_count = collect(counter(inpts))
      ot_count = collect(counter(otpts))
      Reaction(k[t], [S[i] for i in unique(inpts)],
                     [S[o] for o in unique(otpts)],
                     in_count, ot_count)
    end

    ReactionSystem(rxs, t, S, k; name=:model)
end

function simulate(model, u0, p, tspan)
    rxn = MakeReactionSystem(model)
    prob = ODEProblem(rxn, u0, tspan, p)
    sol = solve(prob, Tsit5())
    return sol, prob
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


function generateData(model::AbstractLabelledPetriNet, p, u0, tspan, num_samples, obs_func)
    sample_times = range(tspan[1], stop=tspan[2], length=num_samples)

    sol, prob = simulate(model, u0, p, tspan) # , tstops=sample_times

    obs_true, obs_lbls = obs_func(sol, sample_times)
    
    obs_noisy = obs_true .* (ones(size(obs_true)...)+0.1rand(size(obs_true)...)-0.05ones(size(obs_true)...))

    return obs_noisy, sample_times, prob, sol, obs_true, obs_lbls
end

function makeLossFunction(obs_func, states_to_count, op, tend, sample_data, sample_times, lbls_data)
    # Loss on all populations
    #=loss = function (p)
        sol = solve(remake(op,tspan=(0.0,tend),p=p), Tsit5(), tstops=sample_times)
        vals = hcat(map(ts -> sol.u[findfirst(sol.t .>= ts)], sample_times[1:findlast(sample_times .<= tend)])...)    
        loss = sum(abs2, vals .- sample_vals[:,1:size(vals)[2]])   
        return loss, sol
    end=#
    
    # obs_func = makeObsFunction(proj1_of_pb)

    function loss(p)
        sol = solve(remake(op,tspan=(0.0,tend),p=p), Tsit5(), tstops=sample_times)

        vals, lbls = obs_func(sol, sample_times[1:findlast(sample_times .<= tend)])
        l = 0
        for curr_st in states_to_count
          idx_vals = findfirst(x->x == curr_st,lbls)
          idx_samps = findfirst(x->x == curr_st, lbls_data)
          l += sum(abs2, vals[:,idx_vals] .- sample_data[1:size(vals)[1],idx_samps])
        end

        return l, sol
    end
    return loss
end

function calibrate(ss::StrataSpec, obs_func, states_to_count, u0, p_init, training_data, sample_times, data_labels; name="")
  return calibrate(dom(ss.tpn), obs_func, states_to_count, u0, p_init, training_data, sample_times, data_labels, name=name)
end

function calibrate(model::AbstractPetriNet, obs_func, states_to_count, u0, p_init, training_data, sample_times, data_labels; name="")
    rxn = MakeReactionSystem(model)

    tspan = (sample_times[1], sample_times[end])
    
    prob = ODEProblem(rxn, u0, tspan, p_init)
    
    p_estimate = p_init
    loss = 0
    sol_estimate = Nothing

    for i in 1:(length(sample_times)/10):length(sample_times) # 2:8:50 # 25:25:50 # 
        l_func = makeLossFunction(obs_func, states_to_count, prob, sample_times[Int(floor(i))], training_data, sample_times, data_labels)
        
        res_ode = DiffEqFlux.sciml_train(l_func, p_estimate, #=ADAM(0.05), cb=callback,=# maxiters=2000, abstol=1e-4, reltol=1e-4,
            lower_bounds=repeat([1e-6], length(p_estimate)), upper_bounds=ones(length(p_estimate)))
        p_estimate = res_ode.minimizer
        loss = res_ode.minimum
        
        sol_estimate = solve(remake(prob,tspan=tspan,p=p_estimate), Tsit5())
    end
    
    return  p_estimate, sol_estimate, loss
end

function calibrateFromFiles(fname_mdl,fname_obs_func_and_params,fname_data)
    include(fname_obs_func_and_params)

    mdl = read_json_acset_schema(fname_mdl)

    df_data = CSV.read(fname_data, DataFrame)
    obs_times = 1:length(df_data[:, :t])
    obs_data = df_data
    obs_lbls = names(df_data)

    p_est, sol_est, loss = calibrate(mdl, obs_func, u0, p_init, obs_data, obs_times, obs_lbls; name="")

    # output p_est to file?

    return p_est
end

function plot_obs_w_ests(sample_times, sample_data, sol, obs_func, model::AbstractLabelledPetriNet)
  obs_ests, obs_lbls = obs_func(sol, sample_times)
  plot(sample_times, sample_data, seriestype=:scatter,label="")
  plot!(sample_times,obs_ests, lw=2#=, label=obs_lbls=#)
end

# function stratifyFromFiles() end
# function testCalibrateFromFiles() end
# function testStratifyFromFiles() end
