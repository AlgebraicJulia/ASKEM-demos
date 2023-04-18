# Function to assemble a multi-region model
function make_rgn(n,f_aug=false)
    lstates = []
    ltrans = []
    for ii in 1:n
        push!(lstates,Symbol("Rgn"*string(ii)))
    end
    for ii in 1:n
        for jj in 1:n
            if ii != jj
                push!(ltrans,Symbol("travel"*string(ii)*string(jj)) => ((lstates[ii])=>(lstates[jj])))
            end
        end
    end
    if f_aug
        for ii in 1:n
            push!(ltrans,Symbol("id_d"*string(ii)) => ((lstates[ii])=>(lstates[ii])))
        end
    end
    if f_aug
        for ii in 1:n
            push!(ltrans,Symbol("id_i"*string(ii)) => ((lstates[ii],lstates[ii])=>(lstates[ii],lstates[ii])))
        end
    end
    map(println,ltrans)
    MultiRgn = LabelledPetriNet(lstates,ltrans...)  
end

# Multi-region augmented model
MultiRgn_aug = make_rgn(num_rgns, true)
AlgebraicPetri.Graph(MultiRgn_aug)

# Typed MR model
mr_S = repeat([s],ns(MultiRgn_aug))
mr_T = repeat([t_strata],num_rgns*(num_rgns-1))
mr_I = repeat([i_strata],num_rgns*(num_rgns-1))
mr_O = repeat([o_strata],num_rgns*(num_rgns-1))
for ii in 1:num_rgns
    push!(mr_T,t_disease)
    push!(mr_I,i_disease)
    push!(mr_O,o_disease)
end
for ii in 1:num_rgns
    push!(mr_T,t_interact)
    push!(mr_I,i_interact1)
    push!(mr_I,i_interact2)
    push!(mr_O,o_interact1)
    push!(mr_O,o_interact2)
end
MR_aug_typed = ACSetTransformation(MultiRgn_aug, types,
	  S = mr_S,
	  T = mr_T,
	  I = mr_I,
	  O = mr_O,
	  Name = name -> nothing # specify the mapping for the loose ACSet transform
)
@assert is_natural(MR_aug_typed)
AlgebraicPetri.Graph(dom(MR_aug_typed))

# Function to form stratified model
typed_stratify(typed_model1, typed_model2) =
		Theories.compose(proj1(CategoricalAlgebra.pullback(typed_model1, typed_model2)), typed_model1)

# Multi-region SIRD stratified model
SIRD_MR = typed_stratify(SIRD_aug_typed, MR_aug_typed)
AlgebraicPetri.Graph(dom(SIRD_MR))


function make_living(n,f_aug=false)
    lstates = []
    ltrans = []
    for ii in 1:n
        push!(lstates,Symbol("Live"*string(ii)))
    end
    if f_aug
        for ii in 1:n
            push!(ltrans,Symbol("l_id_d"*string(ii)) => ((lstates[ii])=>(lstates[ii])))
        end
    end
    if f_aug
        for ii in 1:n
            push!(ltrans,Symbol("l_id_s"*string(ii)) => ((lstates[ii])=>(lstates[ii])))
        end
    end
    for ii in 1:n
        for jj in 1:n
            push!(ltrans,Symbol("l_interact"*string(ii)*string(jj)) => ((lstates[ii],lstates[jj])=>(lstates[ii],lstates[jj])))
        end
    end
    map(println,ltrans)
    MultiRgn = LabelledPetriNet(lstates,ltrans...)  
end

# Living augmented model
Living_aug = make_living(num_rgns, true)
AlgebraicPetri.Graph(Living_aug)

# Typed Living model
l_S = repeat([s],ns(Living_aug))
l_T = repeat([t_disease],num_rgns)
l_I = repeat([i_disease],num_rgns)
l_O = repeat([o_disease],num_rgns)
for ii in 1:num_rgns
    push!(l_T,t_strata)
    push!(l_I,i_strata)
    push!(l_O,o_strata)
end
for ii in 1:num_rgns
    for jj in 1:num_rgns
        push!(l_T,t_interact)
        push!(l_I,i_interact1)
        push!(l_I,i_interact2)
        push!(l_O,o_interact1)
        push!(l_O,o_interact2)
    end
end
Living_aug_typed = ACSetTransformation(Living_aug, types,
	  S = l_S,
	  T = l_T,
	  I = l_I,
	  O = l_O,
	  Name = name -> nothing # specify the mapping for the loose ACSet transform
)
@assert is_natural(Living_aug_typed)
AlgebraicPetri.Graph(dom(Living_aug_typed))

simple_trip_typed = typed_stratify(MR_aug_typed, Living_aug_typed)
# AlgebraicPetri.Graph(dom(simple_trip_typed))
   
SIRD_trip = typed_stratify(SIRD_aug_typed,simple_trip_typed) 
# AlgebraicPetri.Graph(dom(SIRD_trip))



function make_nbr_rgn(ladj,f_aug=false)
    n = length(ladj)
    lstates = []
    ltrans = []
    for ii in 1:n
        push!(lstates,Symbol("Rgn"*string(ii)))
    end
    if f_aug
        for ii in 1:n
            push!(ltrans,Symbol("id_s"*string(ii)) => ((lstates[ii])=>(lstates[ii])))
        end
    end
    if f_aug
        for ii in 1:n
            push!(ltrans,Symbol("id_d"*string(ii)) => ((lstates[ii])=>(lstates[ii])))
        end
    end
    for ii in 1:n
        push!(ltrans,Symbol("i_local"*string(ii)) => ((lstates[ii],lstates[ii])=>(lstates[ii],lstates[ii])))
        for jj in ladj[ii]
            push!(ltrans,Symbol("i_travel"*string(ii)*string(jj)) => ((lstates[ii],lstates[jj])=>(lstates[ii],lstates[jj])))
        end
    end
    map(println,ltrans)
    nbr_rgn = LabelledPetriNet(lstates,ltrans...)  
end

# Multi-region augmented model
num_rgns = 4
NbrRgn_aug = make_nbr_rgn([[2,3],[1,4],[1,4],[2,3]], true)
AlgebraicPetri.Graph(NbrRgn_aug)

# Typed MR model
mr_S = repeat([s],num_rgns)
mr_T = repeat([t_strata],num_rgns)
mr_I = repeat([i_strata],num_rgns)
mr_O = repeat([o_strata],num_rgns)
for ii in 1:num_rgns
    push!(mr_T,t_disease)
    push!(mr_I,i_disease)
    push!(mr_O,o_disease)
end
for ii in 1:num_rgns
    for jj in 1:(length(ladj[ii])+1)
        push!(mr_T,t_interact)
        push!(mr_I,i_interact1)
        push!(mr_I,i_interact2)
        push!(mr_O,o_interact1)
        push!(mr_O,o_interact2)
    end
end
Nbr_aug_typed = ACSetTransformation(NbrRgn_aug, types,
	  S = mr_S,
	  T = mr_T,
	  I = mr_I,
	  O = mr_O,
	  Name = name -> nothing # specify the mapping for the loose ACSet transform
)
@assert is_natural(Nbr_aug_typed)
AlgebraicPetri.Graph(dom(Nbr_aug_typed))

# Function to form stratified model
typed_stratify(typed_model1, typed_model2) =
		Theories.compose(proj1(CategoricalAlgebra.pullback(typed_model1, typed_model2)), typed_model1)

# Multi-region SIRD stratified model
SIRD_MR = typed_stratify(SIRD_aug_typed, MR_aug_typed)
AlgebraicPetri.Graph(dom(SIRD_MR))


# *******
using AlgebraicPetri.TypedPetri

const infectious_ontology = LabelledPetriNet(
  [:Pop],
  :infect=>((:Pop, :Pop)=>(:Pop, :Pop)),
  :disease=>(:Pop=>:Pop),
  :strata=>(:Pop=>:Pop)
)

Graph(infectious_ontology)