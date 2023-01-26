
using Catlab, Catlab.CategoricalAlgebra, Catlab.Programs, Catlab.WiringDiagrams, Catlab.Graphics
using AlgebraicPetri
using AlgebraicPetri.BilayerNetworks
using AlgebraicDynamics.UWDDynam
using LabelledArrays
using OrdinaryDiffEq, DelayDiffEq
using Plots

using ASKEM.Dec2022Demo: formSIRD, formInfType, augLabelledPetriNet, sirdAugStates, typeSIRD, 
                      makeMultiAge, typeAge, typed_stratify, formVax, vaxAugStates, typeVax, writeMdlStrat,
                      loadSVIIvR, sviivrAugStates, typeSVIIvR

#*****
# Q2 *
#*****

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

# SEIRD
function formSEIRDnat()
  SEIRDnat = LabelledPetriNet([:S, :E, :I, :R, :D],
  :inf => ((:S, :I)=>(:E, :I)),
  :conv => (:E=>:I),
  :rec => (:I=>:R),
  :death => (:I=>:D),
  :nat_d_s => (:S=>:D),
  :nat_d_e => (:E=>:D),
  :nat_d_i => (:I=>:D),
  :nat_d_r => (:R=>:D),
)
  return SEIRDnat
end

seirdnat = formSEIRDnat()
seirdnat_aug = augLabelledPetriNet(seirdnat,[:S, :E, :I, :R])
seirdnat_typed = ACSetTransformation(seirdnat_aug, types,
  S = [s, s, s, s, s],
  T = [t_interact, t_disease, t_disease, t_disease, t_disease, t_disease, t_disease, t_disease, t_strata, t_strata, t_strata, t_strata],
  I = [i_interact1, i_interact2, i_disease, i_disease, i_disease, i_disease, i_disease, i_disease, i_disease, i_strata, i_strata, i_strata, i_strata],
  O = [o_interact1, o_interact2, o_disease, o_disease, o_disease, o_disease, o_disease, o_disease, o_disease, o_strata, o_strata, o_strata, o_strata],
  Name = name -> nothing 
  )
@assert is_natural(seirdnat_typed)

# Vax
vax_lpn = formVax()
vax_aug_st = vaxAugStates()
vax_aug = augLabelledPetriNet(vax_lpn,vax_aug_st)
vax_typed = typeVax(vax_aug,types)

# Stratified
seirdnat_vax = typed_stratify(seirdnat_typed, vax_typed)

#=function formSEIRHD()
  SEIRHD = LabelledPetriNet([:S, :E, :I, :R, :H, :D],
  :inf => ((:S, :I)=>(:E, :I)),
  :conv => (:E=>:I),
  :rec => (:I=>:R),
  :hosp => (:I=>:H),
  :death => (:H=>:D)
)
  return SEIRHD
end=#

#*****
# Q3 *
#*****

# SEIRDnat "stratified with vax"
function formSEIRDnatV()
  SEIRDnatV = LabelledPetriNet([:Sv, :Ev, :Iv, :Rv, :D],
  :inf => ((:Sv, :Iv)=>(:Ev, :Iv)),
  :conv => (:Ev=>:Iv),
  :rec => (:Iv=>:Rv),
  :death => (:Iv=>:D),
  :nat_d_s => (:Sv=>:D),
  :nat_d_e => (:Ev=>:D),
  :nat_d_i => (:Iv=>:D),
  :nat_d_r => (:Rv=>:D),
)
  return SEIRDnatV
end

SEIRD_composition_pattern = @relation (S, E, I, R, Sv, Ev, Iv, Rv, D) where (S, E, I, R, Sv, Ev, Iv, Rv, D) begin
  SEIRDnat(S, E, I, R, D)
  vax(Sv, Ev, Iv, Rv, D)
  cross_exposure(S, E, I, Sv, Ev, Iv)
end

seirdnat = formSEIRDnat()
seirdnat_v = formSEIRDnatV()

cross_exposure = Open(LabelledPetriNet([:S, :E, :I, :Sv, :Ev, :Iv],
  :inf_uv => ((:S, :Iv) => (:E, :Iv)),
  :inf_vu => ((:Sv, :I) => (:Ev, :I))
))

seirdnat_2x = oapply(SEIRD_composition_pattern, Dict(
  :SEIRDnat => Open(seirdnat),
  :vax => Open(seirdnat_v),
  :cross_exposure => cross_exposure
)) |> apex

# CHIMESVIIvR
sviivr_lbn = read_json_acset(LabelledBilayerNetwork,"../data/CHIME_SVIIvR_dynamics_BiLayer.json")
sviivr = LabelledPetriNet()
migrate!(sviivr,sviivr_lbn)

#***
# Q3b
#***
# Pairwise comparison
max_12 = mca(dom(seirdnat_vax), seirdnat_2x)
max_13 = mca(dom(seirdnat_vax), sviivr)
max_23 = mca(seirdnat_2x, sviivr)

# Three-way comparison

# Plot structural comparisons

#***
# Q3c
#***
# Sim all three models
# i)  Sim w/ vax effic 75%, pop vax 10%
# ii) Sim w/ vax effic 75%, pop vax 10%
# Compare sim outputs

sig(x) = 1/(1+exp(-x))
invsig(x) = log(x/(1-x))

function formTVParams(p,alpha,t_start)
  p_t = [new_p(u,t)=pp for pp in p]
  new_p(u,t) = p[2]*(1 - sig(alpha)*sig(t-t_start))
  p_t[2] = new_p(u,t) 
  new_p(u,t) = p[3]*(1 - sig(alpha)*sig(t-t_start))
  p_t[3] = new_p(u,t) 
  new_p(u,t) = p[4]*(1 - sig(alpha)^2*sig(t-t_start))
  p_t[4] = new_p(u,t) 
  return p_t
end

kappa = 0.1
u0 = [.999*(1-kappa),0,.001,0,0,.999*kappa,0,0,0,0]
p_fixed = [0.000025,0.005,0.001]
tspan = (1,100)
alpha_policy = invsig(0.75) # vax effectiveness
tstart_policy = 1
alpha_init = [0.0]

tv_rxn = vectorfield(dom(seirdnat_vax))
p_t = formTVParams(p,alpha_policy,tstart_policy)
tv_prob = ODEProblem(tv_rxn, u0, tspan, p_t)
tv_sol = solveODE(tv_prob)

kappa = 0.8
u0 = [.999*(1-kappa),0,.001,0,0,.999*kappa,0,0,0,0]
tv_prob = ODEProblem(tv_rxn, u0, tspan, p_t)
tv_sol = solveODE(tv_prob)


#*****
# Q4 *
#*****
# Create eq-wt form_ensemble_model
# Sim for 3ci and 3cii
# Compare w/ indiv run_model_selection

#*****
# Q5 *
#*****
# Sensitivity sensitivity analysis

#*****
# Q6 *
#*****
# Stratify one model (all models) by age.
age_aug = makeMultiAge(16)
age_typed = typeAge(age_aug,types)
seirdnat_vax_age = typed_stratify(seirdnat_vax, age_typed)

#***
# Q6d
#***
# Sim model(s) for...
# i)   High vax rt (+80%) for 65+yo and low vax rt (<15%) for all others
# ii)  High vax rt for all grps
# iii) Repeat i and ii with 20% decrease in contact for school-age children
# iv)  Compare outputs

