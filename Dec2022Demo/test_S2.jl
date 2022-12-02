using Catlab.CategoricalAlgebra
using Catlab.Present, Catlab.Theories
using AlgebraicPetri
using AlgebraicPetri: Graph

include("LibForOct2022Demo.jl");

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
SIRD_aug = LabelledPetriNet([:S, :I, :R, :D],
	  :inf => ((:S, :I)=>(:I, :I)),
	  :rec => (:I=>:R),
	  :death => (:I=>:D),
	  :id => (:S => :S),
	  :id => (:I => :I),
	  :id => (:R => :R)
	)
AlgebraicPetri.Graph(SIRD_aug)

# Typed SIRD model
SIRD_aug_typed = ACSetTransformation(SIRD_aug, types,
S = [s, s, s, s],
T = [t_interact, t_disease, t_disease, t_strata, t_strata, t_strata],
I = [i_interact1, i_interact2, i_disease, i_disease, i_strata, i_strata, i_strata],
O = [o_interact1, o_interact2, o_disease, o_disease, o_strata, o_strata, o_strata],
Name = name -> nothing # specify the mapping for the loose ACSet transform
)
@assert is_natural(SIRD_aug_typed)
AlgebraicPetri.Graph(dom(SIRD_aug_typed))

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

# Multi-region augmented model
MultiAge_aug = makeMultiAge(num_ages, true)
AlgebraicPetri.Graph(MultiAge_aug)

# Typed MA model
ma_S = repeat([s],ns(MultiAge_aug))
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
	  Name = name -> nothing # specify the mapping for the loose ACSet transform
)
@assert is_natural(MA_aug_typed)
AlgebraicPetri.Graph(dom(MA_aug_typed))

# Function to form stratified model
typed_stratify(typed_model1, typed_model2) =
		Theories.compose(proj1(CategoricalAlgebra.pullback(typed_model1, typed_model2)), typed_model1)

# Multi-region SIRD stratified model
SIRD_MA = typed_stratify(SIRD_aug_typed, MA_aug_typed)
AlgebraicPetri.Graph(dom(SIRD_MA))


# Vaccination augmented model
Vax_aug = LabelledPetriNet([:U, :V],
	  :infuu => ((:U, :U)=>(:U, :U)),
	  :infvu => ((:V, :U)=>(:V, :U)),
	  :infuv => ((:U, :V)=>(:U, :V)),
	  :infvv => ((:V, :V)=>(:V, :V)),
	  :vax => (:U => :V),
	  :id => (:U => :U),
	  :id => (:V => :V)
	)
AlgebraicPetri.Graph(Vax_aug)

# Typed Vax model
Vax_aug_typed = ACSetTransformation(Vax_aug, types,
S = [s, s],
T = [t_interact, t_interact, t_interact, t_strata, t_disease, t_disease],
I = [i_interact1, i_interact2, i_interact1, i_interact2, i_interact1, i_interact2, i_strata, i_disease, i_disease],
O = [o_interact1, o_interact2, o_interact1, o_interact2, o_interact1, o_interact2, o_strata, o_disease, o_disease],
Name = name -> nothing # specify the mapping for the loose ACSet transform
)
@assert is_natural(Vax_aug_typed)
AlgebraicPetri.Graph(dom(Vax_aug_typed))

# Vax Multi-region SIRD stratified model
SIRD_MA_Vax = typed_stratify(SIRD_MA, Vax_aug_typed)
AlgebraicPetri.Graph(dom(SIRD_MA_Vax))
