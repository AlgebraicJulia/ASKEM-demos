using ASKEM.Upstream: presentationToLabelledPetriNet, deserialize_wiringdiagram, vectorfield
using ASKEM.Dec2022Demo: formSIRD, formTVParams, solveODE, zeroVal, runControlOptim, makeK, runControlAuto, draw, sig, invsig



function formSEIRHD()
  SEIRHD = LabelledPetriNet([:S, :E, :I, :R, :H, :D],
  :inf => ((:S, :I)=>(:E, :I)),
  :conv => (:E=>:I),
  :rec => (:I=>:R),
  :hosp => (:I=>:H),
  :death => (:H=>:D)
)
  return SEIRHD
end

seirhd = formSEIRHD()

seirhd_typed = ACSetTransformation(seirhd, types,
  S=[s, s, s, s, s, s],
  T=[t_interact, t_disease, t_disease, t_disease, t_disease],
  I=[i_interact1, i_interact2, i_disease, i_disease, i_disease, i_disease],
  O=[o_interact1, o_interact2, o_disease, o_disease, o_disease, o_disease],
  Name=name -> nothing
)
@assert is_natural(seirhd_typed)


vax_lpn = LabelledPetriNet([:U, :V],
  :infuu => ((:U, :U) => (:U, :U)),
  :infvu => ((:V, :U) => (:V, :U)),
  :infuv => ((:U, :V) => (:U, :V)),
  :infvv => ((:V, :V) => (:V, :V)),
  :vax => (:U => :V),
)
# vax_aug = augLabelledPetriNet(vax_lpn, vax_aug_st)

Vax_aug_typed = ACSetTransformation(vax_lpn, types,
  S=[s, s],
  T=[t_interact, t_interact, t_interact, t_interact, t_strata],
  I=[i_interact1, i_interact2, i_interact1, i_interact2, i_interact1, i_interact2, i_interact1, i_interact2, i_strata],
  O=[o_interact1, o_interact2, o_interact1, o_interact2, o_interact1, o_interact2, o_interact1, o_interact2, o_strata],
  Name=name -> nothing
)
@assert is_natural(Vax_aug_typed)

seirhd_vax = stratify_typed(
  seirhd_typed=>[[:strata],[:strata],[:strata],[:strata],[:strata],[]],
  Vax_aug_typed=>[[:disease,:infect],[:disease,:infect]], 
  types′)

@assert is_natural(seirhd_vax)

AlgebraicPetri.Graph(dom(seirhd_vax))



  sig(x) = 1/(1+exp(-x))
  invsig(x) = log(x/(1-x))
  
  #function formTVParams(p,alpha,t_start)
  #function formTVParams(p,alpha,start_time_stop_time)
  function formTVParams(p,alpha,mask_transmission_fraction)
      #t_start = first(start_time_stop_time)
      #t_end = last(start_time_stop_time)
      t_start = tspan[1] + (tspan[2] - tspan[1])*sig(alpha[1])
      t_end = t_start + exp(alpha[2])
      sig_sharpness = 5
      new_p_i(u,t) = p[1] * (1 - mask_transmission_fraction*sig(sig_sharpness*(t-t_start)) + mask_transmission_fraction*(1-sigmoid(sig_sharpness*(t-t_end)))) 
      new_p_r(u,t) = p[2]
      new_p_d(u,t) = p[3]
      p_t = [new_p_i,new_p_r,new_p_d]
      return p_t
  end
  

#function runControlOptim(SIRD, tv_prob, p, tspan, hosp_rt, thresh_H, t_start, alpha_init)
# tspan is the tspan over entire time series
function runControlOptim(SIRD, tv_prob, p, tspan, hosp_rt, thresh_H, t_start, mask_transmission_fraction, alpha_init)
    t = tspan[1]:tspan[2]
    
    function obsHfromI(model::AbstractLabelledPetriNet, sol, sample_times, hosp_rt)
        hosp_sample_vals = hosp_rt * sumvarsbyname(model, :I, sol, sample_times)
        labels = ["H"]
        return hosp_sample_vals, labels
    end
    
    # This is a different alpha.
    # alpha will be the vector of start time and stop times.
    # alpha is the knobs
    # t_start should no longer be an explicit parameter
    function predictH(alpha)
    #function predictH(start_time_stop_time)
        #new_p_t = formTVParams(p,alpha,tspan)
        #new_p_t = formTVParams(p,alpha,start_time_stop_time)
        new_p_t = formTVParams(p,alpha,mask_transmission_fraction)
        tv_sol = solve(remake(tv_prob,tspan=tspan,p=new_p_t), Tsit5(), tstops=t)    
        vals, _ = obsHfromI(SIRD, tv_sol, t[1:findlast(t .<= t[end])], hosp_rt)
        return vals'
    end

    # loss function
    function l_func(control_params)
        # Bias against being above the target hospitalization rate
        pred_H = predictH(control_params)
        breach = sum(pred_H .> thresh_H)/length(pred_H)
        breach_bias = 20
        # Heavily bias against a "negative duration" i.e. end is before start.
        duration_of_masking = exp(control_params[2])# - control_params[1]
        #duration_timespan = tspan[2]
        duration_timespan = tspan[2] - tspan[1]
        #proportion_of_time_masked = duration_of_masking / duration_timespan
        mask_policy_bias = 20
        sig_sharpness = 5
        #new_p_i(u,t) = p[1] * (1 - alpha*sig(sig_sharpness*(t-t_start)) + alpha*(1-sigmoid(sig_sharpness*(t-t_end)))) 
        #tspan[1] + exp(control_params[2])*sig(p[1])t_start
        # Take exp of the duration of masking as a fraction of the total time span
        #l = exp(control_params[2]) + exp(20*breach) - 1
        #l = sig(control_params[1]) + exp(20*breach) - 1
        #l = exp(mask_policy_bias*proportion_of_time_masked) + exp(breach_bias*breach) - 1
        l = (mask_policy_bias*duration_of_masking) + exp(breach_bias*breach) - 1
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

u0 = [999,1,0,0]
p_fixed = [0.000025,0.005,0.001]
tspan = (1,1000)
alpha_policy = invsig(0.05)
tstart_policy = 1
hosp_rt = 0.5
thresh_H = 125
t_start = 200
alpha_init = [0.0]

## Wiring diagram for control as optimization
#s1_sird_cntrl_optim = @program ControlWorkflow (u0::StateVect,p::RateVect,tspan::TSpan,hosp_rt::Rate,thresh_H::Thresh,tstart::TimePt,alpha_init::RateFrac) begin 
SIRD = formSIRD()
tv_rxn = vectorfield(SIRD)
zero_alpha = zeroVal()
#p_t = formTVParams(p,zero_alpha,tstart)
p_t = formTVParams(p,zero_alpha,0.5)
tv_prob = ODEProblem(tv_rxn, u0, tspan, p_t)
alpha, tv_sol_cntrl, t, obs_hosp = runControlOptim( SIRD, tv_prob, p, tspan, hosp_rt, thresh_H, tstart, alpha_init)
return alpha, tv_sol_cntrl, t, obs_hosp
#end



#*****
# Q1 *
#*****
# Forecast cases and hospitalizations for (??) model for three months, no intervention
#   Possibly modify model first
#   Possibly fit model first

#*****
# Q2 *
#*****
# Probability hospitalizations stay under threshold

#*****
# Q3 *
#*****
# Assume consistent social distancing/masking results in 50% decrease in transmission(?). 
# Minimize time the policy is in place. (Can’t be restarted once ended.) 
# From Dec. 28, what are the optimal start and end dates to keep hospitalizations below threshold over the three-month period? 
# How many fewer hospitalizations and cases are there?

#*****
# Q4 *
#*****
# Separate of Q3, sim where policies kick-in when hospitalizations rise above 80% of threshold and stop when fall below.
# Q4a) When are the policies expected to first kick-in?

#***
# Q4b
#***
# Minimize impact on transmission (in terms of % decrease) the policies need to have the first time to...
# i)  Ensure won’t reach threshold at any time in three-month period?
# ii) Ensure the policies only need to be implement once and never reimplemented?

#*****
# Q5 *
#*****
# Instead of NPIs, minimize change in vaccinations needed to have the same impact on cases and hospitalizations as answer from Q3? 
# (Depending on model, may be increase in total vax population or daily vax rate, etc.)

