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

using ASKEM
using Test
using ASKEM.Ontologies: infectious_ontology, strip_names
using ASKEM.Dec2022Demo: formAugSIR, formAugSIRD, formAugSIRD2, formAugQuarantine, altTypeSIR, altTypeSIRD, altTypeSIRD2, altTypeQuarantine, formTarget, formModelList
using ASKEM.Upstream: presentationToLabelledPetriNet, deserialize_wiringdiagram

# Form ComparisonWorkflow presentation of FreeBiproductCategory
@present ComparisonWorkflow(FreeBiproductCategory) begin
    (MdlAug,MdlTyped,MdlTypedList)::Ob 

    formAugSIR::Hom(munit(),MdlAug)
    formAugSIRD::Hom(munit(),MdlAug)
    formAugSIRD2::Hom(munit(),MdlAug)
    formAugQuarantine::Hom(munit(),MdlAug)
    
    altTypeSIR::Hom(MdlAug,MdlTyped)
    altTypeSIRD::Hom(MdlAug,MdlTyped)
    altTypeSIRD2::Hom(MdlAug,MdlTyped)
    altTypeQuarantine::Hom(MdlAug,MdlTyped)

    formTarget::Hom(MdlTyped⊗MdlTyped,MdlTyped)
    formModelList::Hom(MdlTyped⊗MdlTyped⊗MdlTyped⊗MdlTyped,MdlTypedList)

    decompose::Hom(MdlTyped⊗MdlTypedList,MdlTypedList)    
end

# Wiring diagram of comparison example
find_sird_q_components = @program ComparisonWorkflow () begin 
    SIR = formAugSIR()
    SIRD = formAugSIRD()
    SIRD2 = formAugSIRD2()
    Quarantine = formAugQuarantine()
    SIR_typed = altTypeSIR(SIR)
    SIRD_typed = altTypeSIRD(SIRD)
    SIRD2_typed = altTypeSIRD2(SIRD2)
    Quarantine_typed = altTypeQuarantine(Quarantine)

    tgt_mdl = formTarget(SIRD_typed,Quarantine_typed)
    mdl_list = formModelList(SIR_typed,SIRD_typed,Quarantine_typed,SIRD2_typed)

    res = decompose(tgt_mdl,mdl_list)

    return res
end


odirpath = joinpath(@__DIR__,"../outputs")
try mkdir(odirpath) catch end

write_json_acset(find_sird_q_components.diagram,joinpath(odirpath,"s3_find_sird_q.json"))
cwf_lpn = presentationToLabelledPetriNet(ComparisonWorkflow)
write_json_acset(cwf_lpn,joinpath(odirpath,"s3_compar_wf_present.json"))

find_rt = deserialize_wiringdiagram(joinpath(odirpath,"s3_find_sird_q.json"))
find_comps_hom_expr = to_hom_expr(FreeBiproductCategory, find_rt)
find_comps_jfunc = Catlab.Programs.GenerateJuliaPrograms.compile_expr(find_comps_hom_expr)
res = eval(find_comps_jfunc)()
#=
@test res == [(SIRD_typed => Quarantine_typed),
              (Quarantine_typed => SIRD2_typed)]
Footer
=#

AlgebraicPetri.Graph(res[1][1])
AlgebraicPetri.Graph(res[1][2])
AlgebraicPetri.Graph(res[2][1])
AlgebraicPetri.Graph(res[2][2])