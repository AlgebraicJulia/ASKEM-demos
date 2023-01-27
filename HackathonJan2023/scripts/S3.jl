#*****
# Q1 *
#*****
# Forecast cumulative cases and cumulative deaths for MechBayes (SEIRD?) model, no intervention. Provide 90% prediction intervals
# Compare accuracy to true data over span

function formSEIRHD()
    SEIRHD = LabelledPetriNet([:S, :E, :I, :R, :H, :D],
	  :exp => ((:S, :I)=>(:E, :I)),
	  :conv => (:E=>:I),
	  :rec => (:I=>:R),
	  :hosp => (:I=>:H),
      :death => (:H=>:D),
	)
    return SEIRHD
end

seirhd = formSEIRHD()
write_json_acset(seirhd, joinpath("outputs/mdl_jsons/", "s3_seirhd.json"))

function formSEIRD()
    SEIRD = LabelledPetriNet([:S, :E, :I, :R, :H, :D],
	  :exp => ((:S, :I)=>(:E, :I)),
	  :conv => (:E=>:I),
	  :rec => (:I=>:R),
	  :hosp => (:I=>:H),
      :death => (:H=>:D),
	)
    return SEIRD
end

seird = formSEIRD()
write_json_acset(seird, joinpath("outputs/mdl_jsons/", "s3_seird.json"))

function formSIRHD()
    SIRHD = LabelledPetriNet([:S, :E, :I, :R, :H, :D],
	  :exp => ((:S, :I)=>(:E, :I)),
	  :conv => (:E=>:I),
	  :rec => (:I=>:R),
	  :hosp => (:I=>:H),
      :death => (:H=>:D),
	)
    return SIRHD
end

sirhd = formSIRHD()
write_json_acset(sirhd, joinpath("outputs/mdl_jsons/", "s3_sirhd.json"))


valueat(x::Number, u, t) = x
valueat(f::Function, u, t) = try f(u,t) catch e f(t) end

""" vectorfield(pn::AbstractPetriNet)

Generates a Julia function which calculates the vectorfield of the Petri net
being simulated under the law of mass action.

The resulting function has a signature of the form `f!(du, u, p, t)` and can be
passed to the DifferentialEquations.jl solver package.
"""
vectorfield_rw(pn::AbstractPetriNet) = begin
  tm = TransitionMatrices(pn)
  dt = tm.output - tm.input
  f(du,u,p,t) = begin
    rates = zeros(eltype(du),nt(pn)*2)
    u_m = [u[sname(pn, i)] for i in 1:ns(pn)] # [u[i] for i in 1:ns(pn)] # 
    p_m = [u[tname(pn, i)] for i in 1:nt(pn)]
    for i in 1:nt(pn)
      rates[i] = valueat(p_m[i],u,t) * prod(u_m[j] ^ tm.input[i,j] for j in 1:ns(pn))
    end
    #=for i in (nt(pn)+1):(2*nt(pn))
      rates[i] = 0
    end=#
    for j in 1:ns(pn)
      du[sname(pn, j)] = sum(rates[i] * dt[i,j] for i in 1:nt(pn); init = 0.0)
    end
    for j in 1:nt(pn)
      du[tname(pn, j)] = 0
    end
    return du
  end
  return f
end

drift = vectorfield_rw(seirhd)

function dispersion(du, u, p, t)
    du[7,1] = 0.0001
end

u0 = [0.999, 0, 0.001, 0, 0, 0, 1.5, 0.25, 0.5, 0.01, 1/25]
# [0.7, 0.01, 0.001]
tspan = (0.0, 10.0)
p = [0]
prob_sde = SDEProblem(drift, dispersion, u0, tspan, p, noise_rate_prototype=zeros(1,1))
ensembleprob = EnsembleProblem(prob_sde)
data = solve(ensembleprob, SOSRI(); saveat=0.1, trajectories=1000)
plot(EnsembleSummary(data))



sig(x) = 1/(1+exp(-x))
invsig(x) = log(x/(1-x))

function formTVParams(p,alpha,t_start)
    new_p_i(u,t) = p[1] * (1 - sig(alpha)*sig(t-t_start)) 
    new_p_r(u,t) = p[2]
    new_p_d(u,t) = p[3]
    p_t = [new_p_i,new_p_r,new_p_d]
    return p_t
end

function solveODE(tv_prob)
    tv_sol = solve(tv_prob, Tsit5())
end

function zeroVal()
    return 0
end

function runControlOptim(SIRD, tv_prob, p, tspan, hosp_rt, thresh_H, t_start, alpha_init)
    t = tspan[1]:tspan[2]
    
    function obsHfromI(model::AbstractLabelledPetriNet, sol, sample_times, hosp_rt)
        hosp_sample_vals = hosp_rt * sumvarsbyname(model, :I, sol, sample_times)
        labels = ["H"]
        return hosp_sample_vals, labels
    end
    
    function predictH(alpha)
        new_p_t = formTVParams(p,alpha,t_start)
        tv_sol = solve(remake(tv_prob,tspan=tspan,p=new_p_t), Tsit5(), tstops=t)    
        vals, _ = obsHfromI(SIRD, tv_sol, t[1:findlast(t .<= t[end])], hosp_rt)
        return vals'
    end

    function l_func(control_params)
        pred_H = predictH(control_params[1])
        breach = sum(pred_H .> thresh_H)/length(pred_H)
        l = sig(control_params[1]) + exp(20*breach) - 1
        return l
    end

    optf = Optimization.OptimizationFunction((x, p) -> l_func(x), Optimization.AutoForwardDiff())
    optprob = Optimization.OptimizationProblem(optf, alpha_init) 
    result1 = Optimization.solve(optprob, ADAM(0.1),  maxiters = 300) 
    # optprob2 = remake(optprob,u0 = result1.u)
    # result2 = Optimization.solve(optprob2, Optim.BFGS(initial_stepnorm=0.01)) 

    alpha = result1.u[1]
    new_p_t = formTVParams(p,alpha,t_start)
    tv_sol_cntrl = solve(remake(tv_prob,tspan=tspan,p=new_p_t), Tsit5())
    obs_hosp = tv_sol_cntrl(t)[2,:]*hosp_rt

    return alpha, tv_sol_cntrl, t, obs_hosp
end

# Turing
u0 = [1.0, 1.0]
tspan = (0.0, 10.0)
function multiplicative_noise!(du, u, p, t)
    x, y = u
    du[1] = p[5] * x
    return du[2] = p[6] * y
end
p = [1.5, 1.0, 3.0, 1.0, 0.1, 0.1]

function lotka_volterra!(du, u, p, t)
    x, y = u
    α, β, γ, δ = p
    du[1] = dx = α * x - β * x * y
    return du[2] = dy = δ * x * y - γ * y
end

prob_sde = SDEProblem(lotka_volterra!, multiplicative_noise!, u0, tspan, p)

ensembleprob = EnsembleProblem(prob_sde)
data = solve(ensembleprob, SOSRI(); saveat=0.1, trajectories=1000)
plot(EnsembleSummary(data))

@model function fitlv_sde(data, prob)
    # Prior distributions.
    σ ~ InverseGamma(2, 3)
    α ~ truncated(Normal(1.3, 0.5); lower=0.5, upper=2.5)
    β ~ truncated(Normal(1.2, 0.25); lower=0.5, upper=2)
    γ ~ truncated(Normal(3.2, 0.25); lower=2.2, upper=4)
    δ ~ truncated(Normal(1.2, 0.25); lower=0.5, upper=2)
    ϕ1 ~ truncated(Normal(0.12, 0.3); lower=0.05, upper=0.25)
    ϕ2 ~ truncated(Normal(0.12, 0.3); lower=0.05, upper=0.25)

    # Simulate stochastic Lotka-Volterra model.
    p = [α, β, γ, δ, ϕ1, ϕ2]
    predicted = solve(prob, SOSRI(); p=p, saveat=0.1)

    # Early exit if simulation could not be computed successfully.
    if predicted.retcode !== :Success
        Turing.@addlogprob! -Inf
        return nothing
    end

    # Observations.
    for i in 1:length(predicted)
        data[:, i] ~ MvNormal(predicted[i], σ^2 * I)
    end

    return nothing
end;

model_sde = fitlv_sde(odedata, prob_sde)

setadbackend(:forwarddiff)
chain_sde = sample(
    model_sde,
    NUTS(0.25),
    5000;
    init_params=[1.5, 1.3, 1.2, 2.7, 1.2, 0.12, 0.12],
    progress=false,
)
plot(chain_sde)


#*****
# Q2 *
#*****
# What is probability hospitalizations stay under deaths threshold for interval

#*****
# Q3 *
#*****
# If institute mask req, what is prob of staying below threshold. 
#  May need to incorporate masking
#  May want to minimize decrease to meet threshold
#  Use TA1 to translate intervention into transmission parameter, including uncertainty

#*****
# Q4 *
#*****
# Detection rate param, varies as Gauss RW on log-odds scale, likely captures other changes as well
# a) If increase detection by 20% (no other interventions), does that increase forecasted cumulative cases and deaths?
#    How does increase in detection affect the uncertainty of estimates
#    Characterize relationship b/w detection rate, forecasts, uncertainties
#    Does improving detection rates provide decision makers with more accurate forecasts or narrower prediction intervals
# b) Compute accuracy of forecasts (no mask policy) and determine if detection rate improves accuracy

#*****
# Q5 *
#*****
# Model modification, model space exploration, and model run_model_selection
# a) Convert SEIRHD to SEIRHD
#    Compute forecasts and compare accuracy to those of 1
# b) Further modify and do model space exploration and run_model_selection
#    Fit to previous month. Select from forecast accuracy
# c) Three-way structural comparison

#*****
# Q6 *
#*****
# Modify SEIRD to SEIRD+Renewal
# Explain and plot how models differ

#*****
# Q7 *
#*****
# Latest mask mandate can be imposed to ensure with 90% prob won't exceed thresh 
# Characterize relationship between extra cumulative deaths for each day delay