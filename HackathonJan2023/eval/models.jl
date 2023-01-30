#*****
# S1 *
#*****
#= sir_lbn = read_json_acset(LabelledBilayerNetwork,"../data/CHIME_SIR_dynamics_BiLayer.json")
sir = LabelledPetriNet()
migrate!(sir,sir_lbn)
=#

function form_sir()
    SIR = LabelledPetriNet([:S, :I, :R],
    :inf => ((:S, :I)=>(:I, :I)),
    :rec => (:I=>:R),
  )
    return SIR
  end
sir = form_sir()

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


function make_multi_age(n;f_aug=true)
    lstates = []
    ltrans = []
    for ii in 1:n
        push!(lstates,Symbol("Age"*string(ii)))
    end
    for ii in 1:n
        for jj in 1:n
            push!(ltrans,Symbol(("inf",string(ii),string(jj))) => ((lstates[ii],lstates[jj])=>(lstates[ii],lstates[jj])))
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
    return MultiAge
end


age_aug = make_multi_age(3)
age_typed = typeAge(age_aug,types)
sir_age3 = typed_stratify(sir_typed, age_typed)
write_json_acset(dom(sir_age3),"sir_age3.json")
Graph(dom(sir_age3))
age_aug = make_multi_age(16)
age_typed = typeAge(age_aug,types)
sir_age16 = typed_stratify(sir_typed, age_typed)
write_json_acset(dom(sir_age16),"sir_age16.json")

#*****************
# Matrices for S1*
#*****************

using LinearAlgebra

function form_uniform_matrix(dim::Int=3)
  #julia> form_uniform_matrix()
  #3×3 Matrix{Float64}:
  # 0.333333  0.333333  0.333333
  # 0.333333  0.333333  0.333333
  # 0.333333  0.333333  0.333333
  fill(1/dim, (dim,dim))
end

function form_significant_ingroup_matrix(dim::Int=3, ingroup_contact::Float64=2/3)
  #julia> form_significant_ingroup_matrix()
  #3×3 Matrix{Float64}:
  # 0.666667  0.333333  0.333333
  # 0.333333  0.666667  0.333333
  # 0.333333  0.333333  0.666667
  M = fill(1/dim, (dim,dim))
  M[diagind(M)] .= ingroup_contact
  M
end

function form_no_intergroup_contact_matrix(dim::Int=3, ingroup_contact::Float64=1/3)
  #julia> form_no_intergroup_contact_matrix()
  #3×3 Matrix{Float64}:
  # 0.333333  0.0       0.0
  # 0.0       0.333333  0.0
  # 0.0       0.0       0.333333
  M = zeros((dim,dim))
  M[diagind(M)] .= ingroup_contact
  M
end

function form_scaled_uniform_matrix(dim::Int=3, scale::Float64=0.5)
  #julia> form_scaled_uniform_matrix()
  #3×3 Matrix{Float64}:
  # 0.166667  0.166667  0.166667
  # 0.166667  0.166667  0.166667
  # 0.166667  0.166667  0.166667
  scale*fill(1/dim, (dim,dim))
end

function form_skewed_compliance_matrix(dim::Int=3, min_contact::Float64=0.1, max_contact::Float64=0.9)
  #julia> form_skewed_compliance_matrix()
  #3×3 Matrix{Float64}:
  # 0.9       0.333333  0.333333
  # 0.333333  0.5       0.333333
  # 0.333333  0.333333  0.1
  M = fill(1/dim, (dim,dim))
  M[diagind(M)] .= range(max_contact, min_contact, length=dim)
  M
end

#*****
# S2 *
#*****
sidarthe = read_json_acset(LabelledPetriNet,"sidarthe.json")

mca_sidarthe_v = mca(sidarthe,sidarthe_v)
AlgebraicPetri.Graph(mca_sidarthe_v[1])

#*****
# S3 *
#*****
function form_sird()
    SIRD = LabelledPetriNet([:S, :I, :R, :D],
    :inf => ((:S, :I)=>(:I, :I)),
    :rec => (:I=>:R),
    :death => (:I=>:D),
  )
    return SIRD
end
function form_sirh()
    SIRH = LabelledPetriNet([:S, :I, :R, :H],
    :inf => ((:S, :I)=>(:I, :I)),
    :rec => (:I=>:R),
    :hosp => (:I=>:H),
    :hrec => (:H=>:R),
  )
    return SIRH
end
function form_sirhd()
    SIRHD = LabelledPetriNet([:S, :I, :R, :H, :D],
    :inf => ((:S, :I)=>(:I, :I)),
    :rec => (:I=>:R),
    :hosp => (:I=>:H),
    :hrec => (:H=>:R),
    :death => (:H=>:D),
  )
    return SIRHD
end

sir = form_sir()
write_json_acset(sir,"sir.json")
sird = form_sird()
write_json_acset(sird,"sird.json")
sirh = form_sirh()
write_json_acset(sirh,"sirh.json")
sirhd = form_sirhd()
write_json_acset(sirhd,"sirhd.json")

vax_lpn = formVax()
vax_aug_st = vaxAugStates()
vax_aug = augLabelledPetriNet(vax_lpn,vax_aug_st)
vax_typed = typeVax(vax_aug,types)

sirhd_aug = augLabelledPetriNet(sirhd,[:S, :I, :R, :H])
sirhd_typed = ACSetTransformation(sirhd_aug, types,
  S = [s, s, s, s, s],
  T = [t_interact, t_disease, t_disease, t_disease, t_disease, t_strata, t_strata, t_strata, t_strata],
  I = [i_interact1, i_interact2, i_disease, i_disease, i_disease, i_disease, i_strata, i_strata, i_strata, i_strata],
  O = [o_interact1, o_interact2, o_disease, o_disease, o_disease, o_disease, o_strata, o_strata, o_strata, o_strata],
  Name = name -> nothing 
  )
@assert is_natural(sirhd_typed)

sirhd_vax = typed_stratify(sirhd_typed, vax_typed)
write_json_acset(dom(sirhd_vax),"sirhd_vax.json")

n = 16
age_aug = make_multi_age(n)
age_typed = typeAge(age_aug,types)
sirhd_vax_age16 = typed_stratify(sirhd_vax, age_typed)
write_json_acset(dom(sirhd_vax_age16),"sirhd_vax_age16.json")

function form_sirhd_renew()
    SIRHD_renew = LabelledPetriNet([:S, :I, :R, :H, :D],
    :inf => ((:S, :I)=>(:I, :I)),
    :rec => (:I=>:R),
    :hosp => (:I=>:H),
    :hrec => (:H=>:R),
    :death => (:H=>:D),
    :renew => (:R=>:S),
  )
    return SIRHD_renew
end

sirhd_renew = form_sirhd_renew()
write_json_acset(sirhd_renew,"sirhd_renew.json")

sirhd_renew_aug = augLabelledPetriNet(sirhd_renew,[:S, :I, :R, :H])
sirhd_renew_typed = ACSetTransformation(sirhd_renew_aug, types,
  S = [s, s, s, s, s],
  T = [t_interact, t_disease, t_disease, t_disease, t_disease, t_disease, t_strata, t_strata, t_strata, t_strata],
  I = [i_interact1, i_interact2, i_disease, i_disease, i_disease, i_disease, i_disease, i_strata, i_strata, i_strata, i_strata],
  O = [o_interact1, o_interact2, o_disease, o_disease, o_disease, o_disease, o_disease, o_strata, o_strata, o_strata, o_strata],
  Name = name -> nothing 
  )
@assert is_natural(sirhd_renew_typed)

sirhd_renew_vax = typed_stratify(sirhd_renew_typed, vax_typed)
write_json_acset(dom(sirhd_renew_vax),"sirhd_renew_vax.json")
sirhd_renew_vax_age16 = typed_stratify(sirhd_renew_vax, age_typed)
write_json_acset(dom(sirhd_renew_vax_age16),"sirhd_renew_vax_age16.json")

