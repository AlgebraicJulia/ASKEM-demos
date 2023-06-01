using ASKEM

using Catlab.CategoricalAlgebra.CSets

using AlgebraicPetri


##############################
# NOTEBOOK
#bplot(x) = AlgebraicPetri.Graph(x)
#=Suppose we are interested in two processes P₁ and P₂ which we suspect (perhaps thanks to our real-world inutitions) to be related 
and suppose that we have modeled these processes with two models M₁ and M₂. Although it might be difficult to investigate the relationship between 
processes P₁ and P₂ (whose description may be mathematically opaque), it is significantly easier to investigate the relationship between Theories
mathematical models M₁ and M₂. Indeed, whenever M₁ and M₂ are represented as acsets, this investigation reduces to the problem of finding Maximum 
Common Submodel. This problem is a special case of the Maxium Common ACSet which is stated formally below. 

Maximum Common Acset (MCA).
Input: two Acsets a₁ and a₂
Task : find all a with with |a| maximum possible such that there is a monic span of Acset a₁ ← a → a₂.

Finding a maximum common acset is easy in catlab: you can simply call mca(M₁, M₂). Let's demonstrate this. 

For concreteness, let's work with Petri nets, specifically Petri nets which represent epidemiological models. 
One such net may be the SIR model shown below.  
=#
sir    = read_json_acset(LabelledPetriNet, "../data/SIR.json")
AlgebraicPetri.Graph(sir)
#=Another such Petri net might be the SIS model. We've drawn it below. 
=#
sis    = read_json_acset(LabelledPetriNet, "../data/SIS.json")
AlgebraicPetri.Graph(sis)
#=Now if we are interested in determining all of the maximum common submodels of teh SIR and SIS models, we can simply run
=#
maxcommon = mca(sir, sis)
#=In this case there are two maximum common submodels; these are:
=#
AlgebraicPetri.Graph(maxcommon[1])
#=and
=#
AlgebraicPetri.Graph(maxcommon[2])


