using ModelingToolkit
using AlgebraicPetri
using AlgebraicPetri.Epidemiology
using AlgebraicPetri.BilayerNetworks

using Catlab, Catlab.Theories
using Catlab.CategoricalAlgebra
using Catlab.Graphics
using Catlab.Graphics: Graphviz
import Catlab.CategoricalAlgebra: migrate!
using Catlab.WiringDiagrams
using Catlab.Programs
using Catlab.Programs.RelationalPrograms
import Catlab.WiringDiagrams.DirectedWiringDiagrams: WiringDiagramACSet
import Catlab.CategoricalAlgebra.CSets: parse_json_acset

# using JSON

using Random
using OrdinaryDiffEq
using Optimization
using OptimizationOptimisers

import ASKEM.Oct2022Demo: sumvarsbyname, MakeReactionSystem

function formSIRD()
  SIRD = LabelledPetriNet([:S, :I, :R, :D],
                          :inf => ((:S, :I)=>(:I, :I)),
                          :rec => (:I=>:R),
                          :death => (:I=>:D),
                          )
  return SIRD
end

using AlgebraicPetri: vectorfield

valueat(x::Number, u, t) = x
valueat(f::Function, u, t) = try f(u,t) catch e f(t) end
AlgebraicPetri.vectorfield(pn::AbstractPetriNet) = begin
  tm = TransitionMatrices(pn)
  dt = tm.output - tm.input
  f(du,u,p,t) = begin
    rates = zeros(eltype(du),nt(pn))
    # u_m = [u[sname(pn, i)] for i in 1:ns(pn)]
    # p_m = [p[tname(pn, i)] for i in 1:nt(pn)]
    u_m = u
    p_m = p
    for i in 1:nt(pn)
      rates[i] = valueat(p_m[i],u,t) * prod(u_m[j] ^ tm.input[i,j] for j in 1:ns(pn))
    end
    for j in 1:ns(pn)
      # du[sname(pn, j)] = sum(rates[i] * dt[i,j] for i in 1:nt(pn); init = 0.0)
      du[j] = sum(rates[i] * dt[i,j] for i in 1:nt(pn); init = 0.0)
    end
    return du
  end
  return f
end

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


using LinearAlgebra
function makeK(u,p)
  A = [1-p[1]*u[2] -p[1]*u[1] 0 0 2*u[1]u[2]*p[1] 0 0 0;
       p[1]*u[2] 1+p[1]*u[1]-p[2]-p[3] 0 0 0 -2*u[1]*u[2]*p[1]+u[2]*p[2]+u[2]*p[3] 0 0;
       0 p[2] 1 0 0 0 -p[2]*u[2] 0;
       0 p[3] 0 1 0 0 0 -p[3]*u[2];
       0 0 0 0 1 0 0 0;
       0 0 0 0 0 1 0 0;
       0 0 0 0 0 0 1 0;
       0 0 0 0 0 0 0 1]
  
  B = [-u[2]*u[3] 0 0;
       u[2]*u[3] -u[2] -u[2];
       0 u[2] 0;
       0 0 u[2];
       0 0 0;
       0 0 0;
       0 0 0;
       0 0 0]
  
  Q = zeros(8,8)
  Q[1,1] = 1
  Q[2,2] = 1
  Q[4,4] = 1
  
  R = Matrix{Float64}(I,3,3)

  Qf = Q
  
  PN = Q + A'*Qf*A - A'*Qf*B*inv(R+B'*Qf*B)*B'*Qf*A
  return K = inv(R + B'*PN*B)*B'*PN*A
end

function runControlAuto(rxn,u0,p,tspan)
  u_prev = u0
  p_prev = p
  t = tspan[1]:tspan[2]
  t_prev = t[1]
  sol_discr = zeros(length(u_prev),length(t))
  sol_discr[:,1] = u_prev
  for ii in 2:length(t)
    t_curr = t[ii]
    K_curr = makeK(u_prev,p_prev)
    u_aug = [u_prev; [1,1,1,1]]
    p_curr = -K_curr*u_aug
    prob = ODEProblem(rxn, u_prev, (t_prev,t_curr), p_curr)
    sol_curr = solve(prob, Tsit5())
    u_curr = sol_curr(t_curr)
    sol_discr[:,ii] = u_curr
    u_prev = u_curr
    p_prev = p_curr
    t_prev = t_curr
  end
  return sol_discr, t
end

# Form ControlWorkflow presentation of FreeBiproductCategory
@present ControlWorkflow(FreeBiproductCategory) begin
  (LPN,TVRXN,RateVect,RateFrac,TimePt,RateFuncVect,StateVect,TSpan,TVODEProb,TVODESol,
   Rate,Thresh,TimeVect,SolDiscr,RXN)::Ob

  formSIRD::Hom(munit(),LPN)
  vectorfield::Hom(LPN,TVRXN)
  formTVParams::Hom(RateVect⊗RateFrac⊗TimePt,RateFuncVect)
  ODEProblem::Hom(TVRXN⊗StateVect⊗TSpan⊗RateFuncVect,TVODEProb)
  solveODE::Hom(TVODEProb,TVODESol)

  zeroVal::Hom(munit(),RateFrac)
  runControlOptim::Hom(LPN⊗TVODEProb⊗RateVect⊗TSpan⊗Rate⊗Thresh⊗TimePt⊗RateFrac,RateFrac⊗TVODESol⊗TimeVect⊗SolDiscr)

  MakeReactionSystem::Hom(LPN,RXN)
  runControlAuto::Hom(RXN⊗StateVect⊗RateVect⊗TSpan,SolDiscr⊗TimeVect)

end

# Wiring diagram of fixed policy control
s1_sird_cntrl_policy = @program ControlWorkflow (u0::StateVect,p::RateVect,tspan::TSpan,alpha::RateFrac,tstart::TimePt) begin 
  SIRD = formSIRD()
  tv_rxn = vectorfield(SIRD)
  p_t = formTVParams(p,alpha,tstart)
  tv_prob = ODEProblem(tv_rxn, u0, tspan, p_t)
  tv_sol = solveODE(tv_prob)
  return tv_sol
end

# Wiring diagram for control as optimization
s1_sird_cntrl_optim = @program ControlWorkflow (u0::StateVect,p::RateVect,tspan::TSpan,hosp_rt::Rate,thresh_H::Thresh,tstart::TimePt,alpha_init::RateFrac) begin 
  SIRD = formSIRD()
  tv_rxn = vectorfield(SIRD)
  zero_alpha = zeroVal()
  p_t = formTVParams(p,zero_alpha,tstart)
  tv_prob = ODEProblem(tv_rxn, u0, tspan, p_t)
  alpha, tv_sol_cntrl, t, obs_hosp = runControlOptim( SIRD, tv_prob, p, tspan, hosp_rt, thresh_H, tstart, alpha_init)
  return alpha, tv_sol_cntrl, t, obs_hosp
end

# Wiring diagram for automated LQR control
s1_sird_cntrl_auto = @program ControlWorkflow (u0::StateVect,p::RateVect,tspan::TSpan) begin 
  SIRD = formSIRD()
  rxn = MakeReactionSystem(SIRD)
  sol_discr, t = runControlAuto(rxn,u0,p,tspan)
  return sol_discr, t
end

#********************************
# Write wiring diagrams to file *
#********************************
write_json_acset(s1_sird_cntrl_policy.diagram, @__DIR__() * "/../outputs/s1_sird_cntrl_policy.json")
write_json_acset(s1_sird_cntrl_optim.diagram, @__DIR__() * "/../outputs/s1_sird_cntrl_optim.json")
write_json_acset(s1_sird_cntrl_auto.diagram, @__DIR__() * "/../outputs/s1_sird_cntrl_auto.json")

#****************************************
# Test functionality of wiring diagrams *
#****************************************
u0 = [999,1,0,0]
p_fixed = [0.000025,0.005,0.001]
tspan = (1,1000)
alpha_policy = invsig(0.05)
tstart_policy = 1
hosp_rt = 0.5
thresh_H = 125
t_start = 200
alpha_init = [0.0]
# Fixed policy - 5% decrease in infection rate on first day
policy_hom_expr = to_hom_expr(FreeBiproductCategory,s1_sird_cntrl_policy)
policy_jfunc = Catlab.Programs.GenerateJuliaPrograms.compile(policy_hom_expr)
policy_tv_sol = policy_jfunc(u0,p_fixed,tspan,alpha_policy,tstart_policy)
# Control via optimization - keep hospitalizations under 125 w/ decrease starting day 200
optim_hom_expr = to_hom_expr(FreeBiproductCategory,s1_sird_cntrl_optim)
optim_jfunc = Catlab.Programs.GenerateJuliaPrograms.compile(optim_hom_expr)
alpha, optim_tv_sol, t, obs_hosp = optim_jfunc(u0,p_fixed,tspan,hosp_rt,thresh_H,t_start,alpha_init)
# Automated LQR control
auto_hom_expr = to_hom_expr(FreeBiproductCategory,s1_sird_cntrl_auto)
auto_jfunc = Catlab.Programs.GenerateJuliaPrograms.compile(auto_hom_expr)
auto_tv_sol, t_auto = auto_jfunc(u0,p_fixed,tspan)

#*********************************************
# Construct presentation as LabelledPetriNet *
#*********************************************
function presentationToLabelledPetriNet(present)
  lpn = LabelledPetriNet(map(Symbol,generators(present,:Ob)))
  state_idx = Dict()
  for ii in 1:ns(lpn)
    state_idx[sname(lpn,ii)] = ii
  end
  num_obs = length(generators(present,:Ob))
  num_homs = length(generators(present,:Hom))
  for curr_hom in generators(present,:Hom)
    i = add_transition!(lpn, tname=Symbol(curr_hom))
    num_args = length(args(dom(curr_hom)))
    if num_args==1
      add_inputs!(lpn,1,i,state_idx[Symbol(dom(curr_hom))])
    else
      add_inputs!(lpn,num_args,repeat([i],num_args),map(x->state_idx[Symbol(x)], args(dom(curr_hom))))
    end
    num_args = length(args(codom(curr_hom)))
    if num_args==1
      add_outputs!(lpn,1,i,state_idx[Symbol(codom(curr_hom))])
    else
      add_outputs!(lpn,num_args,repeat([i],num_args),map(x->state_idx[Symbol(x)], args(codom(curr_hom))))
    end
  end
  return lpn
end

# Form LabelledPetriNet and write to file
cwf_lpn = presentationToLabelledPetriNet(ControlWorkflow)
write_json_acset(cwf_lpn,"s1_cntrl_wf_present.json")

# lpn_rt = read_json_acset(LabelledPetriNet,"s1_cntrl_wf_present.json")
