using ASKEM
using Catlab, Catlab.CategoricalAlgebra
using AlgebraicPetri
using Test

# Form ComparisonWorkflow presentation of FreeBiproductCategory
@present ComparisonWorkflow(FreeBiproductCategory) begin
    (LPN,TypedMdl,TypedMdlList)::Ob 

    formSIR::Hom(munit(),LPN)
    formSIRD::Hom(munit(),LPN)
    formSIRD2::Hom(munit(),LPN)
    formQuarantine::Hom(munit(),LPN)
    
    altTypeSIR::Hom(LPN,TypedMdl)
    altTypeSIRD::Hom(LPN,TypedMdl)
    altTypeSIRD2::Hom(LPN,TypedMdl)
    altTypeQuarantine::Hom(LPN,TypedMdl)

    formTarget::Hom(TypedMdl⊗TypedMdl,TypedMdl)
    formModelList::Hom(TypedMdl⊗TypedMdl⊗TypedMdl⊗TypedMdl,TypedMdlList)

    decompose::Hom(TypedMdl⊗TypedMdlList)    
end

function formTarget(tm1,tm2)
    return first(legs(pullback(tm1, tm2))) ⋅ tm1
end

function formModelList(tm1,tm2,tm3,tm4)
    return [tm1,tm2,tm3,tm4]
end

# Wiring diagram of comparison example
s3_find_sird_q_components = @program ComparisonWorkflow () begin 
    SIR = formSIR()
    SIRD = formSIRD()
    SIRD2 = formSIRD2()
    Quarantine = formQuarantine()
    SIR_typed = altTypeSIR(SIR,TypedMdl)
    SIRD_typed = altTypeSIRD(SIRD,TypedMdl)
    SIRD2_typed = altTypeSIRD2(SIRD2,TypedMdl)
    Quarantine = altTypeQuarantine(Quarantine,TypedMdl)

    tgt_mdl = formTarget(SIRD_typed,Quarantine_Typed)
    mdl_list = formModelList(SIR_typed,SIRD_typed,Quarantine_typed,SIRD2_typed)

    res = decompose(tgt_mdl,mdl_list)

    return res
end

#=
@test res == [(SIRD_typed => Quarantine_typed),
              (Quarantine_typed => SIRD2_typed)]
Footer
=#

write_json_acset(s3_find_sird_q_components.diagram, "s3_find_sird_q.json")
cwf_lpn = presentationToLabelledPetriNet(ComparisonWorkflow)
write_json_acset(cwf_lpn,"s3_compar_wf_present.json")
