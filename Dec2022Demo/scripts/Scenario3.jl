using Catlab, Catlab.Theories, Catlab.CategoricalAlgebra, Catlab.Graphs, Catlab.CategoricalAlgebra.FinSets
using Catlab.Graphics
using Catlab.CategoricalAlgebra.Diagrams
using Catlab.CategoricalAlgebra.FreeDiagrams
using Catlab.ACSetInterface
import ..FinCats: FreeCatGraph, FinDomFunctor, collect_ob, collect_hom
using Catlab.Programs

using AlgebraicPetri

using PartialFunctions
using MLStyle

#############################
# Some utils                 BEGIN
#############################

concatMap(f,xs) = foldl(vcat, map(f, xs), init =[])

# the smaller acset is the "pattern"
function smallFirst(X::ACSet, Y::ACSet) size(X) ≤ size(Y) ? (X,Y) : (Y,X) end

#given an acset X: C→Set and an object c ∈ C, compute all possible X' ↪ X which are given by removing a single element from the set Xc 
function one_removed_subobs(X::ACSet, c)
  C    = acset_schema(X)
  subs = [filter(j -> j != i, parts(X,c)) for i ∈ parts(X,c)]
  mkSubOb(mysub) = hom(
                    ¬¬Subobject( #need double negation to remove dangling edges
                      X, 
                      NamedTuple( 
                        zip(
                          objects(C), 
                          [o == c ? mysub : parts(X, o) for o ∈ objects(C)]
                          ) 
                      ) 
                    )
                  )
  map(mkSubOb, subs)
end
#############################
# Some utils                 END
#############################


"""
Defintion: let 𝐺: C → 𝐒et be a C-set, we define the _size_ of 𝐺 to be ∑_{c ∈ C} |𝐺c|. 
For example, under this definition, the size of: 
  * a graph G is |GE| + |GV| (num edges + num vertices)
  * a Petri net P is |PT| + |PS| + |PI| + |PO| (num transitions + num species + num input arcs + num output arcs).
"""
function size(X::ACSet) foldl(+, [length(parts(X, oₛ)) for oₛ ∈ objects(acset_schema(X))]) end 


"""Get all monomorphisms from an acset X to an acset Y
"""
monos(X::ACSet, Y::ACSet) = homomorphism(X, strip_names(Y); monic = true, type_components=(Name=x->nothing,),)

"""Ask: "does there exists a mono X ↪ Y ?" 
"""
exists_mono(X::ACSet,Y::ACSet)::Bool = is_homomorphic(X, strip_names(Y); monic = true, type_components=(Name=x->nothing,),)


#############################
# Maximum Common Acset (MCA)    BEGIN
#############################
"""Brute-force implementation of Maximum Common Acset (MCA).
Input: two Acsets a₁ and a₂
Task : find all a with with |a| maximum possible such that there is a monic span of Acset a₁ ← a → a₂. 
"""
function mca(XX::ACSet, YY::ACSet)
  (X,Y) = smallFirst(XX,YY)
  exists_mono(X,Y) ? [X] : mca_help(X, Y)
end 

function mca_help(X::ACSet, Y::ACSet)
  C = acset_schema(X) #X: C → Set
  getsmaller(c) = map(dom, one_removed_subobs(X, c))
  # enumerate all sub-acsets X' ↪ X of the acset X: C → Set obtained by removing one point from Xc for some c ∈ C
  oneRemovedFromX = concatMap(getsmaller, filter(c -> !isempty(parts(X,c)), objects(C)))
  #keep only those X' ↪ X such that there exists mono X' ↪ Y
  Y_subs          = filter(χ -> exists_mono(χ, Y), oneRemovedFromX)
  #either terminate or recurse
  if isempty(Y_subs) 
    concatMap(χ -> mca_help(χ, Y), Y_subs)
  else
    ω = maximum(map(size,Y_subs))
    filter(y -> size(y) == ω, Y_subs)
  end
end
#############################
# Maximum Common Acset (MCA)    END
#############################





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
maxcommon_SIR_and_SIS = mca(sir, sis)
#=In this case there are two maximum common submodels; these are:
=#
AlgebraicPetri.Graph(maxcommon[1])
#=and
=#
AlgebraicPetri.Graph(maxcommon[2])


