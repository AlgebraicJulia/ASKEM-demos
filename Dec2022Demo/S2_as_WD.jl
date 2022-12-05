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

num_ages = 2

# Function to assemble a multi-region model
function makeMultiAge(n,f_aug=false)
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
    map(println,ltrans)
    MultiAge = LabelledPetriNet(lstates,ltrans...)  
end

# Multi-age augmented model
MultiAge_aug = makeMultiAge(num_ages, true)
AlgebraicPetri.Graph(MultiAge_aug)

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

# Vax Multi-region SIRD stratified model
SIRD_MA_Vax = typed_stratify(SIRD_MA, Vax_aug_typed)
AlgebraicPetri.Graph(dom(SIRD_MA_Vax))

# Form Workflow presentation of FreeBiproductCategory
@present Workflow(FreeBiproductCategory) begin
    (File,LRN,LPN,TypedLPN,StrataSpec,LPNss,ObsFunc,ParamVec,StateVec,TSpan,NumTS,SampleData,SampleTimes,ODEProb,ODESol,Labels,Loss)::Ob 
    LoadLRN::Hom(File,LRN)
    
    formSIRD::Hom(munit(),MdlAug)
    formVax::Hom(munit(),MdlAug)
    formInfType::Hom(munit(),MdlType)
    makeMultiAge::Hom(NumStrat⊗Bool,MdlAug)

    typeSIRD::Hom(MdlAug⊗MdlType,MdlTyped)
    typeAge::Hom(MdlAug⊗MdlType,MdlTyped)
    typeVax::Hom(MdlAug⊗MdlType,MdlTyped)
    typed_stratify::Hom(MdlTyped⊗MdlTyped,MdlTyped)

    writeMdlStrat(MdlTyped,File)
 end

# Form wiring diagram of load_stratify_calibrate_control Workflow
stratify_sird_age_vax = @program Workflow (num_ages::NumStrat) begin # 
    # Form models
    mdl_sird = formSIRD()
    mdl_vax = formVax()
    mdl_type = formInfType()
    mdl_age = makeMultiAge(num_ages, true)

    # Specify types of models
    mdl_sird_typed = typeSIRD(mdl_sird,mdl_type)
    mdl_age_typed = typeAge(mdl_age,mdl_type)
    mdl_vax_typed = typeVax(mdl_vax,mdl_type)
    
    # Stratify models
    mdl_sird_age = typed_stratify(mdl_sird_typed, mdl_age_typed)
    mdl_sird_age_vax = typed_stratify(mdl_sird_age, mdl_vax_typed)

    # Write stratified model to file
    write_json_acset(dom(mdl_sird_age_vax), "mdl_sird_age_vax.json")

    return  mdl_sird_age_vax
end

# Display wiring diagram of workflow
draw(load_stratify_calibrate_control)

# Write diagram to file as JSON
write_json_acset(load_stratify_calibrate_control.diagram, "diagram_load_strat_calib_cntrl.json")


