#*****
# S1 *
#*****
sir_lbn = read_json_acset(LabelledBilayerNetwork,"../data/CHIME_SIR_dynamics_BiLayer.json")
sir = LabelledPetriNet()
migrate!(sir,sir_lbn)

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

sir_aug = augLabelledPetriNet(sir,[:S, :I, :R])
sir_typed = ACSetTransformation(sir_aug, types,
  S = [s, s, s],
  T = [t_interact, t_disease, t_strata, t_strata, t_strata],
  I = [i_interact1, i_interact2, i_disease, i_strata, i_strata, i_strata],
  O = [o_interact1, o_interact2, o_disease, o_strata, o_strata, o_strata],
  Name = name -> nothing 
  )
@assert is_natural(sir_typed)

age_aug = makeMultiAge(3)
age_typed = typeAge(age_aug,types)
sir_age3 = typed_stratify(sir_typed, age_typed)

write_json_acset(sir_age3,"sir_age3.json")

#*****
# S2 *
#*****
sidarthe = read_json_acset(LabelledPetriNet,"sidarthe.json")

mca_sidarthe_v = mca(sidarthe,sidarthe_v)
AlgebraicPetri.Graph(mca_sidarthe_v[1])

#*****
# S3 *
#*****

vax_lpn = formVax()
vax_aug_st = vaxAugStates()
vax_aug = augLabelledPetriNet(vax_lpn,vax_aug_st)
vax_typed = typeVax(vax_aug,types)

n = 4
age_aug = makeMultiAge(n)
age_typed = typeAge(age_aug,types)
sir_age3 = typed_stratify(sir_typed, age_typed)
write_json_acset(sir_age3,"sir_age"*string(n)*"".json")
