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

import ASKEM.Upstream: presentationToLabelledPetriNet

# Type system model
types′ = LabelledPetriNet([:Pop],
    :infect=>((:Pop, :Pop)=>(:Pop, :Pop)),
    :disease=>(:Pop=>:Pop),
    :strata=>(:Pop=>:Pop))
types = map(types′, Name=name->nothing)

# Parts of type system for ease of reference
s, = parts(types′, :S)
t_interact, t_disease, t_strata = parts(types′, :T)
i_interact1, i_interact2, i_disease, i_strata = parts(types′, :I)
o_interact1, o_interact2, o_disease, o_strata = parts(types′, :O);

# SIRD augmented model
function formSIRD()
    SIRD_aug = LabelledPetriNet([:S, :I, :R, :D],
	  :inf => ((:S, :I)=>(:I, :I)),
	  :rec => (:I=>:R),
	  :death => (:I=>:D),
	  :id => (:S => :S),
	  :id => (:I => :I),
	  :id => (:R => :R)
	)
    return SIRD_aug
end

function formInfType()
  return types
end

# Typed SIRD model
function typeSIRD(SIRD_aug, types)
    SIRD_aug_typed = ACSetTransformation(SIRD_aug, types,
    S = [s, s, s, s],
    T = [t_interact, t_disease, t_disease, t_strata, t_strata, t_strata],
    I = [i_interact1, i_interact2, i_disease, i_disease, i_strata, i_strata, i_strata],
    O = [o_interact1, o_interact2, o_disease, o_disease, o_strata, o_strata, o_strata],
    Name = name -> nothing 
    )
    @assert is_natural(SIRD_aug_typed)
    return SIRD_aug_typed
end

# Function to assemble a multi-region model
function makeMultiAge(n;f_aug=true)
    lstates = []
    ltrans = []
    for ii in 1:n
        push!(lstates,Symbol("Age"*string(ii)))
    end
    for ii in 1:n
        for jj in 1:n
            push!(ltrans,Symbol("inf"*string(jj)*string(ii)) => ((lstates[ii],lstates[jj])=>(lstates[ii],lstates[jj])))
        end
    end
    if f_aug
        for ii in 1:n
            push!(ltrans,Symbol("dis"*string(ii)) => ((lstates[ii])=>(lstates[ii])))
        end
    end
    if f_aug
        for ii in 1:n
            push!(ltrans,Symbol("id"*string(ii)) => ((lstates[ii])=>(lstates[ii])))
        end
    end
    MultiAge = LabelledPetriNet(lstates,ltrans...)
end

# Typed Multi-age model
function typeAge(MultiAge_aug,types)
    num_ages = ns(MultiAge_aug)
    ma_S = repeat([s],num_ages)
    ma_T = repeat([t_interact],num_ages*num_ages)
    ma_I = repeat([i_interact1, i_interact2],num_ages*num_ages)
    ma_O = repeat([o_interact1, o_interact2],num_ages*num_ages)
    for ii in 1:num_ages
        push!(ma_T,t_disease)
        push!(ma_I,i_disease)
        push!(ma_O,o_disease)
    end
    for ii in 1:num_ages
        push!(ma_T,t_strata)
        push!(ma_I,i_strata)
        push!(ma_O,o_strata)
    end
    MA_aug_typed = ACSetTransformation(MultiAge_aug, types,
	  S = ma_S,
	  T = ma_T,
	  I = ma_I,
	  O = ma_O,
	  Name = name -> nothing 
    )
    @assert is_natural(MA_aug_typed)
    return MA_aug_typed
end

# Function to form stratified model
typed_stratify(typed_model1, typed_model2) =
		Theories.compose(proj1(CategoricalAlgebra.pullback(typed_model1, typed_model2)), typed_model1)


# Vaccination augmented model
function formVax()
    Vax_aug = LabelledPetriNet([:U, :V],
	  :infuu => ((:U, :U)=>(:U, :U)),
	  :infvu => ((:V, :U)=>(:V, :U)),
	  :infuv => ((:U, :V)=>(:U, :V)),
	  :infvv => ((:V, :V)=>(:V, :V)),
	  :vax => (:U => :V),
	  :id => (:U => :U),
	  :id => (:V => :V)
	)
    return Vax_aug
end 

# Typed Vax model
function typeVax(Vax_aug,types)
    Vax_aug_typed = ACSetTransformation(Vax_aug, types,
    S = [s, s],
    T = [t_interact, t_interact, t_interact, t_interact, t_strata, t_disease, t_disease],
    I = [i_interact1, i_interact2, i_interact1, i_interact2, i_interact1, i_interact2, i_interact1, i_interact2, i_strata, i_disease, i_disease],
    O = [o_interact1, o_interact2, o_interact1, o_interact2, o_interact1, o_interact2, o_interact1, o_interact2, o_strata, o_disease, o_disease],
    Name = name -> nothing 
    )
    @assert is_natural(Vax_aug_typed)
    return Vax_aug_typed
end

function writeMdlStrat(mdl, file)
  write_json_acset(dom(mdl), file)
end

# Vax Multi-region SIRD stratified model
# SIRD_MA_Vax = typed_stratify(SIRD_MA, Vax_aug_typed)
# AlgebraicPetri.Graph(dom(SIRD_MA_Vax))

# Form Workflow presentation of FreeBiproductCategory
@present StratificationWorkflow(FreeBiproductCategory) begin
    (File,LRN,LPN,TypedLPN,StrataSpec,LPNss,ObsFunc,ParamVec,StateVec,TSpan,NumTS,SampleData,SampleTimes,ODEProb,ODESol,Labels,Loss,MdlAug,MdlType,MdlTyped,NumStrat)::Ob
    LoadLRN::Hom(File,LRN)
    
    formSIRD::Hom(munit(),MdlAug)
    formVax::Hom(munit(),MdlAug)
    formInfType::Hom(munit(),MdlType)
    makeMultiAge::Hom(NumStrat,MdlAug)

    typeSIRD::Hom(MdlAug⊗MdlType,MdlTyped)
    typeAge::Hom(MdlAug⊗MdlType,MdlTyped)
    typeVax::Hom(MdlAug⊗MdlType,MdlTyped)
    typed_stratify::Hom(MdlTyped⊗MdlTyped,MdlTyped)

    writeMdlStrat::Hom(MdlTyped⊗File,munit())
 end

# Form wiring diagram of load_stratify_calibrate_control Workflow
stratify_sird_age_vax = @program StratificationWorkflow (num_ages::NumStrat, out_file::File) begin #
    # Form models
    mdl_sird = formSIRD()
    mdl_vax = formVax()
    mdl_type = formInfType()
    mdl_age = makeMultiAge(num_ages)

    # Specify types of models
    mdl_sird_typed = typeSIRD(mdl_sird,mdl_type)
    mdl_age_typed = typeAge(mdl_age,mdl_type)
    mdl_vax_typed = typeVax(mdl_vax,mdl_type)
    
    # Stratify models
    mdl_sird_age = typed_stratify(mdl_sird_typed, mdl_age_typed)
    mdl_sird_age_vax = typed_stratify(mdl_sird_age, mdl_vax_typed)

    # Write stratified model to file
    writeMdlStrat(mdl_sird_age_vax, out_file)
end

# Display wiring diagram of workflow
# draw(stratify_sird_age_vax)

#********************************
# Write diagram to file as JSON *
#********************************
write_json_acset(stratify_sird_age_vax.diagram, "s2_strat_sird_age_vax.json")


#****************************************
# Test functionality of wiring diagrams *
#****************************************
stratify_sird_hom_expr = to_hom_expr(FreeBiproductCategory, stratify_sird_age_vax)
stratify_sird_jfunc = Catlab.Programs.GenerateJuliaPrograms.compile(stratify_sird_hom_expr)
stratify_sird_jfunc(7,"sird_age7_vax.json")

#******************************
# Confirm read-in of diagrams *
#******************************
import Catlab.WiringDiagrams.DirectedWiringDiagrams: WiringDiagramACSet
using Catlab.CategoricalAlgebra.CSets: parse_json_acset
function Catlab.CategoricalAlgebra.CSets.parse_json_acset(::Type{T}, input::AbstractDict) where T <: ACSet 
    out = T()
    for (k,v) ∈ input
      add_parts!(out, Symbol(k), length(v))
    end
    for l ∈ values(input)
      for (i, j) ∈ enumerate(l)
        for k ∈ keys(j) # (k,v) = j # 
          v = j[k]
          vtype = eltype(out[Symbol(k)])
          if !(v isa vtype)
            v = vtype(v)
          end
          set_subpart!(out, i, Symbol(k), v)
        end
      end
    end
    out
end
  

deserialize_wiringdiagram(filepath::String, value=nothing) = deserialize_wiringdiagram!(read_json_acset(WiringDiagramACSet{Symbol,Any,Any,Any}, filepath), isnothing(value) ? filepath : value)
deserialize_wiringdiagram!(dwd, value) = begin
  convsymbol(dwd, key) = begin
    dwd[key] .= Symbol.(dwd[key])
  end
  
  dwd[:box_type] .= Box{Symbol}
  convsymbol(dwd, :in_port_type)
  convsymbol(dwd, :out_port_type)
  convsymbol(dwd, :outer_in_port_type)
  convsymbol(dwd, :outer_out_port_type)
  dwd[:value] .= map(Symbol,dwd[:value])
  wd_acset2 = WiringDiagramACSet{Symbol,Any,Any,DataType}()
  copy_parts!(wd_acset2,dwd)
  return WiringDiagram{ThBiproductCategory, Symbol, Any, Any}(wd_acset2, value)
end
rt_wd = deserialize_wiringdiagram("s2_strat_sird_age_vax.json")


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
write_json_acset(swf_lpn,"s2_strat_wf_present.json")

# lpn_rt = read_json_acset(LabelledPetriNet,"s1_cntrl_wf_present.json")
