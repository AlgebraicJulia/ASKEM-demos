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
using ASKEM.Dec2022Demo: formTVParams, solveODE, zeroVal, runControlOptim, makeK, runControlAuto, draw, sig, invsig

using Catlab.Present
using AlgebraicPetri: Graph

using ASKEM
using Test
using ASKEM.Ontologies: infectious_ontology, strip_names
using ASKEM.Dec2022Demo: formAugSIR, formAugSIRD, formAugSIRD2, formAugQuarantine, altTypeSIR, altTypeSIRD, altTypeSIRD2, altTypeQuarantine, formTarget, formModelList

using ASKEM.Dec2022Demo: formSIRD, formInfType, augLabelledPetriNet, sirdAugStates, typeSIRD, 
                      makeMultiAge, typeAge, typed_stratify, formVax, vaxAugStates, typeVax, writeMdlStrat,
                      loadSVIIvR, sviivrAugStates, typeSVIIvR



#***********************
# OVERALL PRESENTATION *
#***********************
@present UberWorkflow(FreeBiproductCategory) begin
    (File,BLN,LPN,LPNType,AugStates,LPNAug,LPNTyped,ListLPN,NumStrat,
    RXN,TVRXN,RateVect,RateFuncVect,StateVect,TSpan,TVODEProb,TVODESol,
    SpaceLPN,ObsFunc,LossFunc,SelecFunc,TrainData,TestData,CntrlArgVect,CntrlParamVect,Loss,TimeVect,SolDiscr)::Ob 
    
    #*****
    # S1 *
    #*****
    load_bln::Hom(File,BLN)
    load_tpn::Hom(File,LPNTyped)
    convert_bln_to_lpn::Hom(BLN,LPN)
    
    form_aug_states_s1m1::Hom(munit(),AugStates)
    form_aug_states_s1m2::Hom(munit(),AugStates)
    form_aug_states_s1m3::Hom(munit(),AugStates)
    
    augment_lpn::Hom(LPN⊗AugStates,LPNAug)
    
    form_type_inf::Hom(munit(),LPNType)
    
    type_s1m1::Hom(LPNAug⊗LPNType,LPNTyped)
    type_s1m2::Hom(LPNAug⊗LPNType,LPNTyped)
    type_s1m3::Hom(LPNAug⊗LPNType,LPNTyped)

    stratify_typed_lpn::Hom(LPNTyped⊗LPNTyped,LPNTyped)

    form_vax::Hom(munit(),LPN)
    form_aug_states_vax::Hom(munit(),AugStates)
    type_vax::Hom(LPNAug⊗LPNType,LPNTyped)
    
    form_age::Hom(NumStrat,LPNAug)
    type_age::Hom(LPNAug⊗LPNType,LPNTyped)
    
    mca::Hom(ListLPN,ListLPN)
    # plot_mca::Hom(ListLPN⊗ListLPN)

    MakeReactionSystem::Hom(LPN,RXN)
    vectorfield::Hom(LPN,TVRXN)
    form_tv_params_const::Hom(RateVect,RateFuncVect)
    ODEProblem::Hom(TVRXN⊗StateVect⊗TSpan⊗RateFuncVect,TVODEProb)
    solveODE::Hom(TVODEProb,TVODESol)
    
    fit_model::Hom(LPN⊗ObsFunc⊗LossFunc⊗StateVect⊗RateVect⊗TrainData⊗TimeVect,RateVect⊗TVODESol⊗Loss)

    # form_ensemble_model
    # sensitivity_analysis
    # compare_trajectories

    #*****
    # S2 *
    #*****
    modify_lpn::Hom(LPN,LPN)
    modify_typed_lpn::Hom(LPNTyped,LPNTyped)

    form_tv_params_s2q3::Hom(RateVect⊗CntrlParamVect,RateFuncVect)
    form_tv_params_s2q4::Hom(RateVect⊗CntrlParamVect,RateFuncVect)
    form_tv_params_s2q5::Hom(RateVect⊗CntrlParamVect,RateFuncVect)
    
    run_control_optim_s2q3::Hom(LPN⊗TVODEProb⊗RateVect⊗TSpan⊗CntrlArgVect⊗CntrlParamVect,CntrlParamVect⊗TVODESol⊗TimeVect⊗SolDiscr)  
    run_control_optim_s2q4bi::Hom(LPN⊗TVODEProb⊗RateVect⊗TSpan⊗CntrlArgVect⊗CntrlParamVect,CntrlParamVect⊗TVODESol⊗TimeVect⊗SolDiscr)
    run_control_optim_s2q4bii::Hom(LPN⊗TVODEProb⊗RateVect⊗TSpan⊗CntrlArgVect⊗CntrlParamVect,CntrlParamVect⊗TVODESol⊗TimeVect⊗SolDiscr)
    run_control_optim_s2q5::Hom(LPN⊗TVODEProb⊗RateVect⊗TSpan⊗CntrlArgVect⊗CntrlParamVect,CntrlParamVect⊗TVODESol⊗TimeVect⊗SolDiscr)

    #*****
    # S3 *
    #*****
    form_aug_states_s3m1::Hom(munit(),AugStates)
    form_aug_states_s3m2::Hom(munit(),AugStates)
    form_aug_states_s3m3::Hom(munit(),AugStates)
    form_aug_states_s3m4::Hom(munit(),AugStates)
    
    type_s3m1::Hom(LPNAug⊗LPNType,LPNTyped)
    type_s3m2::Hom(LPNAug⊗LPNType,LPNTyped)
    type_s3m3::Hom(LPNAug⊗LPNType,LPNTyped)
    type_s3m4::Hom(LPNAug⊗LPNType,LPNTyped)
    
    form_mask::Hom(munit(),LPN)
    form_aug_states_mask::Hom(munit(),AugStates)
    type_mask::Hom(LPNAug⊗LPNType,LPNTyped)

    run_model_selection::Hom(SpaceLPN⊗ObsFunc⊗LossFunc⊗SelecFunc⊗StateVect⊗RateVect⊗TrainData⊗TimeVect,RateVect⊗TVODESol⊗Loss)

    # detection_sensitivity
    
    # zeroVal::Hom(munit(),RateFrac)

    run_control_optim_s3q3::Hom(LPN⊗TVODEProb⊗RateVect⊗TSpan⊗CntrlArgVect⊗CntrlParamVect,CntrlParamVect⊗TVODESol⊗TimeVect⊗SolDiscr)
    run_control_optim_s3q7::Hom(LPN⊗TVODEProb⊗RateVect⊗TSpan⊗CntrlArgVect⊗CntrlParamVect,CntrlParamVect⊗TVODESol⊗TimeVect⊗SolDiscr)
    
    #*****
    # S4 *
    #*****
    form_aug_states_s4m1::Hom(munit(),AugStates)
    type_s4m1::Hom(LPNAug⊗LPNType,LPNTyped)

    form_test::Hom(munit(),LPN)
    form_aug_states_test::Hom(munit(),AugStates)
    type_test::Hom(LPNAug⊗LPNType,LPNTyped)
    
    form_cohorts::Hom(munit(),LPN)
    form_aug_states_cohorts::Hom(munit(),AugStates)
    type_cohorts::Hom(LPNAug⊗LPNType,LPNTyped)

    form_tv_params_s4q3a::Hom(RateVect⊗CntrlParamVect,RateFuncVect)
    form_tv_params_s4q4::Hom(RateVect⊗CntrlParamVect,RateFuncVect)

    run_control_optim_s4q3a::Hom(LPN⊗TVODEProb⊗RateVect⊗TSpan⊗CntrlArgVect⊗CntrlParamVect,CntrlParamVect⊗TVODESol⊗TimeVect⊗SolDiscr)
    run_control_optim_s4q4::Hom(LPN⊗TVODEProb⊗RateVect⊗TSpan⊗CntrlArgVect⊗CntrlParamVect,CntrlParamVect⊗TVODESol⊗TimeVect⊗SolDiscr)

    #*****
    # SC *
    #*****
    restratify_typed_lpn::Hom(LPNTyped⊗LPNTyped,LPNTyped)
    
end


import Catlab.Graphics.Graphviz: run_graphviz
draw_diagram(d,s::String;fmt="svg") = open(s, "w") do fp
    run_graphviz(fp, to_graphviz(d),format=fmt)
end

#***********************
# Scenario 1 Workflows *
#***********************

odirpath_fmt = joinpath(@__DIR__,"outputs/wd_figures")
fmt = ".svg"
try mkdir(odirpath_fmt) catch end

#***********************
# Scenario 2 Workflows *
#***********************
s2_q1 = @program UberWorkflow (fname::File,u0::StateVect,p::RateVect,tspan::TSpan) begin 
    s2m1_bln = load_bln(fname)
    s2m1 = convert_bln_to_lpn(s2m1_bln)
    tv_rxn = vectorfield(s2m1)
    p_t = form_tv_params_const(p)
    tv_prob = ODEProblem(tv_rxn, u0, tspan, p_t)
    tv_sol = solveODE(tv_prob)
    return tv_sol
end

s2_q3 = @program UberWorkflow (fname::File,u0::StateVect,p::RateVect,tspan::TSpan,beta::CntrlArgVect,alpha_init::CntrlParamVect) begin 
    s2m3_bln = load_bln(fname)
    s2m3 = convert_bln_to_lpn(s2m3_bln)
    tv_rxn = vectorfield(s2m3)
    p_t = form_tv_params_s2q3(p,alpha_init)
    tv_prob = ODEProblem(tv_rxn, u0, tspan, p_t)
    alpha, tv_sol_cntrl, t, sol_discr = run_control_optim_s2q3(s2m3, tv_prob, p_t, tspan, beta, alpha_init)  
    return alpha, tv_sol_cntrl, t
end

s2_q4 = @program UberWorkflow (fname::File,u0::StateVect,p::RateVect,tspan::TSpan,alpha_init::CntrlParamVect) begin 
    s2m4_bln = load_bln(fname)
    s2m4 = convert_bln_to_lpn(s2m4_bln)
    tv_rxn = vectorfield(s2m4)
    p_t = form_tv_params_s2q4(p,alpha_init)
    tv_prob = ODEProblem(tv_rxn, u0, tspan, p_t)
    tv_sol = solveODE(tv_prob)
    return tv_sol
end

s2_q4bi = @program UberWorkflow (fname::File,u0::StateVect,p::RateVect,tspan::TSpan,beta::CntrlArgVect,alpha_init::CntrlParamVect) begin 
    s2m4_bln = load_bln(fname)
    s2m4 = convert_bln_to_lpn(s2m4_bln)
    tv_rxn = vectorfield(s2m4)
    p_t = form_tv_params_s2q4(p,alpha_init)
    tv_prob = ODEProblem(tv_rxn, u0, tspan, p_t)
    alpha, tv_sol_cntrl, t, sol_discr = run_control_optim_s2q4bi(s2m4, tv_prob, p_t, tspan, beta, alpha_init)
    return alpha, tv_sol_cntrl, t
end

s2_q4bii = @program UberWorkflow (fname::File,u0::StateVect,p::RateVect,tspan::TSpan,beta::CntrlArgVect,alpha_init::CntrlParamVect) begin 
    s2m4_bln = load_bln(fname)
    s2m4 = convert_bln_to_lpn(s2m4_bln)
    tv_rxn = vectorfield(s2m4)
    p_t = form_tv_params_s2q4(p,alpha_init)
    tv_prob = ODEProblem(tv_rxn, u0, tspan, p_t)
    alpha, tv_sol_cntrl, t, sol_discr = run_control_optim_s2q4bii(s2m4, tv_prob, p_t, tspan, beta, alpha_init)
    return alpha, tv_sol_cntrl, t
end

s2_q5 = @program UberWorkflow (fname::File,u0::StateVect,p::RateVect,tspan::TSpan,beta::CntrlArgVect,alpha_init::CntrlParamVect) begin 
    s2m5_bln = load_bln(fname)
    s2m5 = convert_bln_to_lpn(s2m5_bln)
    tv_rxn = vectorfield(s2m5)
    p_t = form_tv_params_s2q5(p,alpha_init)
    tv_prob = ODEProblem(tv_rxn, u0, tspan, p_t)
    alpha, tv_sol_cntrl, t, sol_discr = run_control_optim_s2q5(s2m5, tv_prob, p_t, tspan, beta, alpha_init)
    return alpha, tv_sol_cntrl, t
end

draw_diagram(s2_q1, joinpath(odirpath_fmt, "s2_q1"*fmt))
draw_diagram(s2_q3, joinpath(odirpath_fmt, "s2_q3"*fmt))
draw_diagram(s2_q4, joinpath(odirpath_fmt, "s2_q4"*fmt))
draw_diagram(s2_q4bi, joinpath(odirpath_fmt, "s2_q4bi"*fmt))
draw_diagram(s2_q4bii, joinpath(odirpath_fmt, "s2_q4bii"*fmt))
draw_diagram(s2_q5, joinpath(odirpath_fmt, "s2_q5"*fmt))

#***********************
# Scenario 3 Workflows *
#***********************

#=
s3_q1 = @program ControlWorkflow (u0::StateVect,p::RateVect,tspan::TSpan,hosp_rt::Rate,thresh_H::Thresh,tstart::TimePt,alpha_init::RateFrac) begin 
    s3m1_bln = load_bln(File)
    s3m1 = convert_bln_to_lpn(s3m1_bln)
    tv_rxn = vectorfield(s3m1)
    zero_alpha = zeroVal()
    p_t = formTVParams(p,zero_alpha,tstart)
    tv_prob = ODEProblem(tv_rxn, u0, tspan, p_t)
    alpha, tv_sol_cntrl, t, obs_hosp = runControlOptim( SIRD, tv_prob, p, tspan, hosp_rt, thresh_H, tstart, alpha_init)
    return alpha, tv_sol_cntrl, t, obs_hosp
end
=#

#***********************
# Scenario 4 Workflows *
#***********************
s4_q3a = @program UberWorkflow (fname::File,u0::StateVect,p::RateVect,tspan::TSpan,beta::CntrlArgVect,alpha_init::CntrlParamVect) begin 
    s4m1_bln = load_bln(fname)
    s4m1 = convert_bln_to_lpn(s4m1_bln)

    mtype = form_type_inf()
 
    aug_states = form_aug_states_s4m1()
    s4m1_aug = augment_lpn(s4m1,aug_states)
    s4m1_typed = type_s4m1(s4m1_aug,mtype)
    
    mtests = form_test()
    tests_aug_states = form_aug_states_test()
    mtests_aug = augment_lpn(mtests,tests_aug_states)
    mtests_typed = type_test(mtests_aug,mtype)
    
    mcohorts = form_cohorts()
    cohorts_aug_states = form_aug_states_cohorts()
    mcohorts_aug = augment_lpn(mcohorts,cohorts_aug_states)
    mcohorts_typed = type_cohorts(mcohorts_aug,mtype)

    s4m1_strat1 = stratify_typed_lpn(s4m1_typed,mcohorts_typed)
    s4m1_strat2 = stratify_typed_lpn(s4m1_strat1,mtests_typed)

    tv_rxn = vectorfield(s4m1_strat2)
    p_t = form_tv_params_s4q3a(p,alpha_init)
    tv_prob = ODEProblem(tv_rxn, u0, tspan, p_t)
    alpha, tv_sol_cntrl, t, sol_discr = run_control_optim_s4q3a(s4m1_strat2, tv_prob, p_t, tspan, beta, alpha_init)
    return alpha, tv_sol_cntrl, t
end


s4_q4 = @program UberWorkflow (fname::File,u0::StateVect,p::RateVect,tspan::TSpan,beta::CntrlArgVect,alpha_init::CntrlParamVect) begin 
    s4m1_bln = load_bln(fname)
    s4m1 = convert_bln_to_lpn(s4m1_bln)

    mtype = form_type_inf()
 
    aug_states = form_aug_states_s4m1()
    s4m1_aug = augment_lpn(s4m1,aug_states)
    s4m1_typed = type_s4m1(s4m1_aug,mtype)
    
    mtests = form_test()
    tests_aug_states = form_aug_states_test()
    mtests_aug = augment_lpn(mtests,tests_aug_states)
    mtests_typed = type_test(mtests_aug,mtype)
    
    mcohorts = form_cohorts()
    cohorts_aug_states = form_aug_states_cohorts()
    mcohorts_aug = augment_lpn(mcohorts,cohorts_aug_states)
    mcohorts_typed = type_cohorts(mcohorts_aug,mtype)

    s4m1_strat1 = stratify_typed_lpn(s4m1_typed,mcohorts_typed)
    s4m1_strat2 = stratify_typed_lpn(s4m1_strat1,mtests_typed)

    tv_rxn = vectorfield(s4m1_strat2)
    p_t = form_tv_params_s4q4(p,alpha_init)
    tv_prob = ODEProblem(tv_rxn, u0, tspan, p_t)
    alpha, tv_sol_cntrl, t, sol_discr = run_control_optim_s4q4(s4m1_strat2, tv_prob, p_t, tspan, beta, alpha_init)
    return alpha, tv_sol_cntrl, t
end

draw_diagram(s4_q3a, joinpath(odirpath_fmt, "s4_q3a"*fmt))
draw_diagram(s4_q4, joinpath(odirpath_fmt, "s4_q4"*fmt))

#*********************
# Challenge Scenario *
#*********************
sc = @program UberWorkflow (fname1::File,fname2::File,u0::StateVect,p::RateVect,tspan::TSpan) begin 
    scm_orig = load_tpn(fname1)    
    new_strat = load_tpn(fname2)    
    scm_restrat = restratify_typed_lpn(scm_orig,new_strat)
    tv_rxn = vectorfield(scm_restrat)
    p_t = form_tv_params_const(p)
    tv_prob = ODEProblem(tv_rxn, u0, tspan, p_t)
    tv_sol = solveODE(tv_prob)
    return tv_sol
end

draw_diagram(sc, joinpath(odirpath_fmt, "sc"*fmt))
