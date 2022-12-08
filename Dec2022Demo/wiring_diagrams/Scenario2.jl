using Catlab, Catlab.Theories
using Catlab.CategoricalAlgebra
using Catlab.Graphics
using Catlab.Graphics: Graphviz
import Catlab.CategoricalAlgebra: migrate!
using Catlab.WiringDiagrams
using Catlab.Programs
using Catlab.Programs.RelationalPrograms

using Catlab.Present
using AlgebraicPetri
using AlgebraicPetri: Graph

using ASKEM.Dec2022Demo: formSIRD, formInfType, augLabelledPetriNet, sirdAugStates, typeSIRD, 
                      makeMultiAge, typeAge, typed_stratify, formVax, vaxAugStates, typeVax, writeMdlStrat,
                      loadSVIIvR, sviivrAugStates, typeSVIIvR
using ASKEM.Upstream: presentationToLabelledPetriNet, deserialize_wiringdiagram


# Vax Multi-region SIRD stratified model
# SIRD_MA_Vax = typed_stratify(SIRD_MA, Vax_aug_typed)
# AlgebraicPetri.Graph(dom(SIRD_MA_Vax))

# Form Workflow presentation of FreeBiproductCategory
@present StratificationWorkflow(FreeBiproductCategory) begin
    (LPN,MdlAug,MdlType,MdlTyped,NumStrat,AugStates,File)::Ob
    
    formSIRD::Hom(munit(),LPN)
    formVax::Hom(munit(),LPN)
    formInfType::Hom(munit(),MdlType)
    makeMultiAge::Hom(NumStrat,MdlAug)

    sirdAugStates::Hom(munit(),AugStates)
    vaxAugStates::Hom(munit(),AugStates)
    augLabelledPetriNet::Hom(LPN⊗AugStates,MdlAug)

    typeSIRD::Hom(MdlAug⊗MdlType,MdlTyped)
    typeAge::Hom(MdlAug⊗MdlType,MdlTyped)
    typeVax::Hom(MdlAug⊗MdlType,MdlTyped)
    typed_stratify::Hom(MdlTyped⊗MdlTyped,MdlTyped)

    writeMdlStrat::Hom(MdlTyped⊗File,munit())

    loadSVIIvR::Hom(File,LPN)
    sviivrAugStates::Hom(munit(),AugStates)
    typeSVIIvR::Hom(MdlAug⊗MdlType,MdlTyped) 
    
    # loadBucky::Hom(File,LPN)
 end

# Form wiring diagram of load_stratify_calibrate_control Workflow
stratify_sird_age_vax = @program StratificationWorkflow (num_ages::NumStrat, out_file::File) begin #
    # Form models
    mdl_sird = formSIRD()
    mdl_vax = formVax()
    mdl_type = formInfType()
    mdl_age_aug = makeMultiAge(num_ages)

    # Augment models
    sird_aug_states = sirdAugStates()
    mdl_sird_aug = augLabelledPetriNet(mdl_sird,sird_aug_states)
    vax_aug_states = vaxAugStates()
    mdl_vax_aug = augLabelledPetriNet(mdl_vax,vax_aug_states)

    # Specify types of models
    mdl_sird_typed = typeSIRD(mdl_sird_aug,mdl_type)
    mdl_age_typed = typeAge(mdl_age_aug,mdl_type)
    mdl_vax_typed = typeVax(mdl_vax_aug,mdl_type)
    
    # Stratify models
    mdl_sird_age = typed_stratify(mdl_sird_typed, mdl_age_typed)
    mdl_sird_age_vax = typed_stratify(mdl_sird_age, mdl_vax_typed)

    # Write stratified model to file
    writeMdlStrat(mdl_sird_age_vax, out_file)

    # return mdl_sird_age_vax
end

# Display wiring diagram of workflow
# draw(stratify_sird_age_vax)

#********************************
# Write diagram to file as JSON *
#********************************
odirpath = joinpath(@__DIR__,"../outputs")
write_json_acset(stratify_sird_age_vax.diagram, joinpath(odirpath,"s2_strat_sird_age_vax.json"))


#****************************************
# Test functionality of wiring diagrams *
#****************************************
stratify_sird_hom_expr = to_hom_expr(FreeBiproductCategory, stratify_sird_age_vax)
stratify_sird_jfunc = Catlab.Programs.GenerateJuliaPrograms.compile(stratify_sird_hom_expr)
test = stratify_sird_jfunc(7,joinpath(odirpath,"sird_age7_vax.json"))

#******************************
# Confirm read-in of diagrams *
#******************************
rt_wd = deserialize_wiringdiagram(joinpath(odirpath,"s2_strat_sird_age_vax.json"))

to_graphviz(rt_wd, labels=true)
to_graphviz(stratify_sird_age_vax, labels=true)

using Test 
# Check equality of read-in wd-acset to original
# we aren't getting exact equality
@testset "Round Trip of WiringDiagram" begin
  @test rt_wd != stratify_sird_age_vax
  @test !is_isomorphic(rt_wd, stratify_sird_age_vax)

  @test is_isomorphic(rt_wd.diagram, stratify_sird_age_vax.diagram)
  @test boxes(rt_wd) == boxes(stratify_sird_age_vax)
  @test wires(rt_wd) == wires(stratify_sird_age_vax)
  for i in parts(rt_wd.diagram, :InPort)
    @test in_wires(rt_wd,5) == in_wires(stratify_sird_age_vax,5)
    @test out_wires(rt_wd,5) == out_wires(stratify_sird_age_vax,5)
  end
end
# Form roundtrip wiring diagram from read-in wd acset
# s2_strat = WiringDiagram{ThBiproductCategory,Any,Any,Any}(rt_wd_acset2,nothing)

#*********************************************
# Construct presentation as LabelledPetriNet *
#*********************************************

# Form LabelledPetriNet and write to file
swf_lpn = presentationToLabelledPetriNet(StratificationWorkflow)
write_json_acset(swf_lpn,joinpath(odirpath,"s2_strat_wf_present.json"))

# lpn_rt = read_json_acset(LabelledPetriNet,"s1_cntrl_wf_present.json")


#=
#**************
stratify_sviivr_age = @program StratificationWorkflow (in_file::File, num_ages::NumStrat, out_file::File) begin #
    # Form models
    mdl_sviivr = loadSVIIvR(in_file)
    mdl_type = formInfType()
    mdl_age_aug = makeMultiAge(num_ages)

    # Augment models
    sviivr_aug_states = sviivrAugStates()
    mdl_sviivr_aug = augLabelledPetriNet(mdl_sviivr,sviivr_aug_states)
    
    # Specify types of models
    mdl_sviivr_typed = typeSVIIvR(mdl_sviivr_aug,mdl_type)
    mdl_age_typed = typeAge(mdl_age_aug,mdl_type)
    
    # Stratify models
    mdl_sviivr_age = typed_stratify(mdl_sviivr_typed, mdl_age_typed)

    # Write stratified model to file
    writeMdlStrat(mdl_sviivr_age, out_file)
end

write_json_acset(stratify_sviivr_age.diagram, joinpath(odirpath,"s2_strat_sviivr_age.json"))

rt_wd = deserialize_wiringdiagram(joinpath(odirpath,"s2_strat_sviivr_age.json"))
stratify_sviivr_hom_expr = to_hom_expr(FreeBiproductCategory, rt_wd)
stratify_sviivr_jfunc = Catlab.Programs.GenerateJuliaPrograms.compile(stratify_sviivr_hom_expr)
stratify_sviivr_jfunc("../data/CHIME_SVIIvR_dynamics_BiLayer.json",7,joinpath(odirpath,"sviivr_age7.json"))
=#