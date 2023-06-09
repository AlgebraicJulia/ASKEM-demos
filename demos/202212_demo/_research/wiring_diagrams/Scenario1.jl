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

# using Random
using OrdinaryDiffEq
using Optimization
using OptimizationOptimisers

import ASKEM.Oct2022Demo: sumvarsbyname, MakeReactionSystem

using ASKEM.Upstream: presentationToLabelledPetriNet, deserialize_wiringdiagram, vectorfield
using ASKEM.Dec2022Demo: formSIRD, formTVParams, solveODE, zeroVal, runControlOptim, makeK, runControlAuto, draw, sig, invsig

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
import Catlab.Graphics.Graphviz: run_graphviz
draw_diagram(d,s::String) = open(s, "w") do fp
    run_graphviz(fp, to_graphviz(d))
end

odirpath_dot = joinpath(@__DIR__,"../outputs/s1_diagrams")
try mkdir(odirpath_dot) catch end
draw_diagram(s1_sird_cntrl_policy, joinpath(odirpath_dot, "sird_cntrl_policy.dot"))
draw_diagram(s1_sird_cntrl_optim, joinpath(odirpath_dot, "sird_cntrl_optim.dot"))
draw_diagram(s1_sird_cntrl_auto, joinpath(odirpath_dot, "sird_cntrl_auto.dot"))

odirpath_wd = joinpath(@__DIR__,"../outputs")
try mkdir(odirpath_wd) catch end
write_json_acset(s1_sird_cntrl_policy.diagram, joinpath(odirpath_wd, "s1_sird_cntrl_policy.json"))
write_json_acset(s1_sird_cntrl_optim.diagram, joinpath(odirpath_wd, "s1_sird_cntrl_optim.json"))
write_json_acset(s1_sird_cntrl_auto.diagram, joinpath(odirpath_wd,"s1_sird_cntrl_auto.json"))

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

# Form LabelledPetriNet and write to file
cwf_lpn = presentationToLabelledPetriNet(ControlWorkflow)
write_json_acset(cwf_lpn,joinpath(odirpath_wd,"s1_cntrl_wf_present.json"))

# lpn_rt = read_json_acset(LabelledPetriNet,"s1_cntrl_wf_present.json")
