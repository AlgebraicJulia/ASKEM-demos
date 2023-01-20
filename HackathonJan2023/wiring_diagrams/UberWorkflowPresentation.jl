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


@present UberWorkflow(FreeBiproductCategory) begin
    (File,BLN,LPN,LPNType,AugStates,LPNAug,LPNTyped,ListLPN,
    RXN,TVRXN,RateVect,RateFuncVect,StateVect,TSpan,TVODEProb,TVODESol,
    SpaceLPN,ObsFunc,LossFunc,SelecFunc,TrainData,TestData,CntrlArgVect,CntrlParamVect,Loss,TimeVect,SolDiscr)::Ob 
    
    #*****
    # S1 *
    #*****
    load_bln::Hom(File,BLN)
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
    plot_mca::Hom(ListLPN⊗ListLPN)

    MakeReactionSystem::Hom(LPN,RXN)
    vectorfield::Hom(LPN,TVRXN)
    form_tv_params_const::Hom(RateVect,RateFuncVect)
    ODEProblem::Hom(TVRXN⊗StateVect⊗TSpan⊗RateFuncVect,TVODEProb)
    solveODE::Hom(TVODEProb,TVODESol)
    
    fit_model::Hom(LPN⊗ObsFunc⊗LossFunc⊗StateVect⊗RateVect⊗TrainData⊗TimeVect,RateVect⊗ODESol⊗Loss)

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

    run_model_selection::Hom(SpaceLPN⊗ObsFunc⊗LossFunc⊗SelecFunc⊗StateVect⊗RateVect⊗TrainData⊗TimeVect,RateVect⊗ODESol⊗Loss)

    # detection_sensitivity
    
    zeroVal::Hom(munit(),RateFrac)

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

    run_control_optim_s4q3a::Hom(LPN⊗TVODEProb⊗RateVect⊗TSpan⊗CntrlArgVect⊗CntrlParamVect,CntrlParamVect⊗TVODESol⊗TimeVect⊗SolDiscr)
    run_control_optim_s4q4::Hom(LPN⊗TVODEProb⊗RateVect⊗TSpan⊗CntrlArgVect⊗CntrlParamVect,CntrlParamVect⊗TVODESol⊗TimeVect⊗SolDiscr)

    #*****
    # SC *
    #*****
    restratify_typed_lpn::Hom(LPNTyped,LPNTyped)
    
end
