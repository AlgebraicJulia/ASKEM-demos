#*****
# Q1 *
#*****
# Define test strategy, minimize total testing whil maintaining infections below isolation capacity
# Can have unique test type and frequency per cohort. 
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

# SEIR
function form_seir()
  SEIR = LabelledPetriNet([:S, :E, :I, :R],
  :inf => ((:S, :I)=>(:E, :I)),
  :conv => (:E=>:I),
  :rec => (:I=>:R),
)
  return SEIR
end

#=seiir = formSEIIR()
seiir_aug = augLabelledPetriNet(seiir,[:S, :E, :I, :R])
seiir_typed = ACSetTransformation(seiir_aug, types,
  S = [s, s, s, s, s],
  T = [t_interact, t_disease, t_disease, t_disease, t_disease, t_strata, t_strata, t_strata, t_strata],
  I = [i_interact1, i_interact2, i_disease, i_disease, i_disease, i_disease, i_strata, i_strata, i_strata, i_strata],
  O = [o_interact1, o_interact2, o_disease, o_disease, o_disease, o_disease, o_strata, o_strata, o_strata, o_strata],
  Name = name -> nothing 
  )
@assert is_natural(seiir_typed)
=#

function form_test_iso(st, ed)
    states = [Symbol(st), Symbol("Iso_"*st)]
    if st != ed
        states = push!(states,Symbol(ed))
    end
    # test_and_iso = LabelledPetriNet([Symbol(st), Symbol("Iso_"*st)],
    test_and_iso = LabelledPetriNet(states,
    Symbol("t_rapid_pos_"*st) => (Symbol(st)=>Symbol("Iso_"*st)),
    Symbol("t_rapid_neg_"*st) => (Symbol(st)=>Symbol(st)),
    Symbol("t_pcr_pos_"*st) => (Symbol(st)=>Symbol("Iso_"*st)),
    Symbol("t_pcr_neg_"*st) => (Symbol(st)=>Symbol(st)),
    Symbol("deiso_"*st) => (Symbol("Iso_"*st)=>Symbol(ed))
  )
    return test_and_iso
end


testing_composition_pattern = @relation (S, E, I, R, Iso_s, Iso_e, Iso_i, Iso_r) where (S, E, I, R, Iso_s, Iso_e, Iso_i, Iso_r) begin
  SEIR(S, E, I, R)
  TI_S(S, Iso_s, S)
  TI_E(E, Iso_e, R)
  TI_I(I, Iso_i, R)
  TI_R(R, Iso_r, R)
  # cross_exposure(S, E, I, Sv, Ev, Iv)
end

seir = form_seir()
test_iso_s = form_test_iso("S", "S")
test_iso_e = form_test_iso("E", "R")
test_iso_i = form_test_iso("I", "R")
test_iso_r = form_test_iso("R", "R")

#= cross_exposure = Open(LabelledPetriNet([:S, :E, :I, :Sv, :Ev, :Iv],
  :inf_uv => ((:S, :Iv) => (:E, :Iv)),
  :inf_vu => ((:Sv, :I) => (:Ev, :I))
))=#

seir_test_iso = oapply(testing_composition_pattern, Dict(
  :SEIR => Open(seir),
  :TI_S => Open(test_iso_s),
  :TI_E => Open(test_iso_e),
  :TI_I => Open(test_iso_i),
  :TI_R => Open(test_iso_r)
  # :cross_exposure => cross_exposure
)) |> apex

import Catlab.Graphics.Graphviz: run_graphviz
draw_diagram(d,s::String;fmt="svg") = open(s, "w") do fp
    run_graphviz(fp, to_graphviz(d),format=fmt)
end
draw_diagram(seir_test_iso, joinpath("../outputs/pn_figures/", "seir_test_iso.svg"))

function form_cohort()
    cohort = LabelledPetriNet([:U, :G, :F],
    :inf_uu => ((:U, :U)=>(:U, :U)),
    :inf_ug => ((:U, :G)=>(:U, :G)),
    :inf_uf => ((:U, :F)=>(:U, :F)),
    :inf_gu => ((:G, :U)=>(:G, :U)),
    :inf_gg => ((:G, :G)=>(:G, :G)),
    :inf_gf => ((:G, :F)=>(:G, :F)),
    :inf_fu => ((:F, :U)=>(:F, :U)),
    :inf_fg => ((:F, :G)=>(:F, :G)),
    :inf_ff => ((:F, :F)=>(:F, :F)),
    :dis_u => (:U => :U),
    :dis_g => (:G => :G),
    :dis_f => (:F => :F),
    :strat_u => (:U => :U),
    :strat_g => (:G => :G),
    :strat_f => (:F => :F)
  )
    return cohort
  end

cohort = form_cohort()
cohort_typed = ACSetTransformation(cohort, types,
    S = [s, s, s],
    T = [t_interact, t_interact, t_interact, t_interact, t_interact, t_interact, t_interact, t_interact, t_interact, t_disease, t_disease, t_disease, t_strata, t_strata, t_strata],
    I = [i_interact1, i_interact2, i_interact1, i_interact2, i_interact1, i_interact2, i_interact1, i_interact2, i_interact1, i_interact2, i_interact1, i_interact2, i_interact1, i_interact2, i_interact1, i_interact2, i_interact1, i_interact2, i_disease, i_disease, i_disease, i_strata, i_strata, i_strata],
    O = [o_interact1, o_interact2, o_interact1, o_interact2, o_interact1, o_interact2, o_interact1, o_interact2, o_interact1, o_interact2, o_interact1, o_interact2, o_interact1, o_interact2, o_interact1, o_interact2, o_interact1, o_interact2, o_disease, o_disease, o_disease, o_strata, o_strata, o_strata],
    Name = name -> nothing 
    )
  @assert is_natural(cohort_typed)
  

seir_test_iso_typed = ACSetTransformation(seir_test_iso, types,
  S = [s, s, s, s, s, s, s, s],
  T = [t_interact, t_disease, t_disease, t_strata, t_strata, t_strata, t_strata, t_strata, t_strata, t_strata, t_strata, t_strata, t_strata, t_strata, t_strata, t_strata, t_strata, t_strata, t_strata, t_strata, t_strata, t_strata, t_strata],
  I = [i_interact1, i_interact2, i_disease, i_disease, i_strata, i_strata, i_strata, i_strata, i_strata, i_strata, i_strata, i_strata, i_strata, i_strata, i_strata, i_strata, i_strata, i_strata, i_strata, i_strata, i_strata, i_strata, i_strata, i_strata],
  O = [o_interact1, o_interact2, o_disease, o_disease, o_strata, o_strata, o_strata, o_strata, o_strata, o_strata, o_strata, o_strata, o_strata, o_strata, o_strata, o_strata, o_strata, o_strata, o_strata, o_strata, o_strata, o_strata, o_strata, o_strata],
  Name = name -> nothing 
  )
@assert is_natural(seir_test_iso_typed)

seir_test_iso_cohort = typed_stratify(seir_test_iso_typed, cohort_typed)

AlgebraicPetri.Graph(dom(seir_test_iso_cohort))

write_json_acset(dom(seir_test_iso_cohort), joinpath("outputs/mdl_jsons/", "s4_seir_test_iso_cohort.json"))

#=
function form_test_strata()
    test_strata = LabelledPetriNet([:rapid, :pcr],
    :inf => ((:S, :I)=>(:E, :I)),
    :conv => (:E=>:I),
    :rec => (:I=>:R),
    :death => (:I=>:D),
  )
    return test_strata
end
=#

#*****
# Q3 *
#*****
# Stratify base model (??) by test type and cohorts
# Fit model and optimize

#*****
# Q4 *
#*****
# Optimize strategy wrt costs (antigen tests 1/5 cost of PCR but ~1/2 as sensitive)