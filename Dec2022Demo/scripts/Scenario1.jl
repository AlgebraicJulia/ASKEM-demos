using ASKEM
using AlgebraicPetri
using Catlab
using Catlab.CategoricalAlgebra
using Plots
using ModelingToolkit
using OrdinaryDiffEq
using Optimization
using Flux
using Optim

types′ = infectious_ontology
types = strip_names(infectious_ontology)

p_i(u,t) = 0.000025 # -0.000025*cos(t/100*2*pi)
# sigmoid(x) = 1/(1-exp(-x))
# p_i(u,t) = 0.00005*(1/(1+exp(-(t-30))))
p_r(u,t) = 0.005+0.005*cos(t/100*4*pi+pi/4)
p_d(u,t) = 0.001+0.001*cos(t/100*8*pi+pi/2)

p_t = [p_i,p_r,p_d]

tv_rxn = vectorfield(SIRD)
tspan = (0,1000)
u0 = [999,1,0,0]
tv_prob = ODEProblem(tv_rxn, u0, tspan, p_t)
tv_sol = solve(tv_prob, Tsit5())
plot(tv_sol)

#= p_i(u,t) = 0.000025-0.000025*cos(t/100*2*pi)
p_r(u,t) = 0.01
p_d(u,t) = 0.002
p_t = [p_i,p_r,p_d]=#


#************************************
# Set up optimization for "control" *
#************************************
# t = range(tspan[1],tspan[2])
t = tspan[1]:tspan[2]
inf_rt = 0.000025
hosp_rt = 0.5
thresh_H = 125
t_start = 200
pinit = [0.0]

function obsHfromI(model::AbstractLabelledPetriNet, sol, sample_times, hosp_rt)
  hosp_sample_vals = hosp_rt * sumvarsbyname(model, :I, sol, sample_times)
  labels = ["H"]
  return hosp_sample_vals, labels
end

sig(x) = 1/(1+exp(-x))
function predictH(alpha)
  new_p_i(u,t) = inf_rt * (1 - sig(alpha)*sig(t-t_start)) # t < t_start ? inf_rt : sig(alpha) * inf_rt # 
  new_p_t = [new_p_i,p_r,p_d]
  tv_sol = solve(remake(tv_prob,tspan=tspan,p=new_p_t), Tsit5(), tstops=t)    
  vals, _ = obsHfromI(SIRD, tv_sol, t[1:findlast(t .<= t[end])], hosp_rt)
  return vals'
end

function l_func(control_params)
  pred_H = predictH(control_params[1])
  breach = sum(pred_H .> thresh_H)/length(pred_H)
  l = sig(control_params[1]) + exp(20*breach) - 1
  # println(breach)
  # println(l)
  return l
end
optf = Optimization.OptimizationFunction((x, p) -> l_func(x), Optimization.AutoForwardDiff())
optprob = Optimization.OptimizationProblem(optf, pinit) # , lcons=[0], ucons=[1]
result1 = Optimization.solve(optprob, ADAM(0.1),  maxiters = 300) # callback = callback,
optprob2 = remake(optprob,u0 = result1.u)
result2 = Optimization.solve(optprob2, Optim.BFGS(initial_stepnorm=0.01)) # , callback=callback
    

new_p_i(u,t) = inf_rt * (1 - sig(result1.u[1])*sig(t-t_start)) # t < t_start ? inf_rt : result1.u[1] * inf_rt # 1/(1+exp(-(t-t_start))
new_p_t = [new_p_i,p_r,p_d]
tv_sol_cntrl = solve(remake(tv_prob,tspan=tspan,p=new_p_t), Tsit5())
plot(t,tv_sol_cntrl(t)')
plot!(t,tv_sol_cntrl(t)[2,:]*.5)


#***********************
# Automated controller *
#***********************
# control = -K*u

p_t = [p_i,p_r,p_d]
tv_rxn = vectorfield(SIRD)
tspan = (0,1000)
u0 = [999,1,0,0]
tv_prob = ODEProblem(tv_rxn, u0, tspan, p_t)
tv_sol = solve(tv_prob, Tsit5())
plot(tv_sol)


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

rxn = MakeReactionSystem(SIRD)
function mk_sol_discr()
  tspan = (0,1000)
  u_prev = [999,1,0,0]
  p_prev = [0.000025, 0.005, 0.001]
  prev_t = t[1]
  sol_discr = zeros(length(u_prev),length(t))
  sol_discr[:,1] = u_prev
  for ii in 2:length(t)
    curr_t = t[ii]
    K_curr = makeK(u_prev,p_prev)
    u_aug = [u_prev; [1,1,1,1]]
    p_curr = -K_curr*u_aug
    prob = ODEProblem(rxn, u_prev, (prev_t,curr_t), p_curr)
    sol_curr = solve(prob, Tsit5())
    u_curr = sol_curr(curr_t)
    sol_discr[:,ii] = u_curr
    u_prev = u_curr
    p_prev = p_curr
    prev_t = curr_t
  end
end


#****

function evalparam(β::Function,t)
  β(t)
end
function evalparam(β::Number,t)
  β
end

# model = SIRD

# rxn = MakeReactionSystem(model)

# rxn = vectorfield(model) # f!(du, u, p, t)

# tspan = (sample_times[1], sample_times[end])

# prob = ODEProblem(rxn, u0, tspan, p_init)

# p_estimate = p_init
# loss = 0
# sol_estimate = Nothing

# AtoB = LabelledPetriNet([:A, :B], :recover => (:A=>:B));
# p_r(t) = 0.05+0.01*cos(t/100*4+pi/4)

# function p(t)
#   [p_r(t)]
# end

p_t = Dict()
p_t[:inf] = p_i
p_t[:recover] = p_r
p_t[:death] = p_d

@present TheoryOrigMIRANet <: SchLabelledReactionNet begin
    MID::AttrType
    MCTX::AttrType
    Template::AttrType
    mira_ids::Attr(S, MID)
    mira_context::Attr(S, MCTX)
    mira_initial_value::Attr(S, Rate) # Concentration) #
    template_type::Attr(T, Template)
    parameter_name::Attr(T, Name)
    parameter_value::Attr(T, Rate)
end
@abstract_acset_type AbstractOrigMIRANet <: AbstractLabelledReactionNet
@acset_type OrigMIRANet(TheoryOrigMIRANet) <: AbstractOrigMIRANet

# *****
# lbn_sir = read_json_acset(LabelledBilayerNetwork,"../../data/CHIME_SIR_dynamics_BiLayer.json")
# lbn_sir = read_json_acset(LabelledBilayerNetwork,"../../data/CHIME_SVIIvR_dynamics_BiLayer.json")
# lbn_sir = read_json_acset(LabelledBilayerNetwork,"../../data/Vucky_dynamics_BiLayer.json")
# lpn_sir = LabelledPetriNet()
# migrate!(lpn_sir,lbn_sir)
# AlgebraicPetri.Graph(lpn_sir)

# itype = read_json_acset(LabelledPetriNet,"../data/infectious_type.json")
