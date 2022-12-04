module Oct2022Demo

export stratify_with_obs_function, simulate, calibrate, MakeReactionSystem, generate_data, 
       get_infected_states, get_susceptible_states,
       get_recovered_states, optimise_p, full_train, obs_SIRD, 
       makeLossFunction, plot_obs_w_ests, full_train_rand

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

using ..Stratify # StrataSpec


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
  
  state_samples = sol(sample_times) # WARNING THIS VARIABLE IS UNUSED
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


function stratify_with_obs_function(ss1, ss2, type_system)
  pb = stratify_span(ss1,ss2,type_system)
  f = first(legs(pb)) # proj1 not defined on this type
  return apex(pb), makeObsFunction(f)
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

function plot_obs_w_ests(sample_times, sample_data, sol, obs_func)
  obs_ests, obs_lbls = obs_func(sol, sample_times)
  plot(sample_times, sample_data, seriestype=:scatter,label="")
  plot!(sample_times,obs_ests, lw=2, label=reshape(map(String,obs_lbls),(1, length(obs_lbls))))
end

 
function stateidx(model, name)
  filter(parts(model, :S)) do i
    model[i, :sname] == name
  end
end
  
function stateidx_stratified(model, name, dim=1)
  filter(parts(model, :S)) do i
    model[i, :sname][dim] == name
  end
end
  
function sumvarsbyname(model, name, sol, sample_times)
  idxs = statidx(model, :I)
  sample_vals = sum(sol(sample_times)[idxs,:], dims=1)
end

function obs_IHD_from_mdl(model::AbstractLabelledPetriNet, sol, sample_times)
  inf_sample_vals = sumvarsbyname(model, :I, sol, sample_times)
  hosp_sample_vals = sumvarsbyname(model, :H, sol, sample_times)
  dead_sample_vals = sumvarsbyname(model, :D, sol, sample_times)

  labels = reshape(["I", "H", "D"],1,3)

  return hcat(inf_sample_vals, hosp_sample_vals, dead_sample_vals), labels
end

# function stratifyFromFiles() end
# function testCalibrateFromFiles() end
# function testStratifyFromFiles() end
end