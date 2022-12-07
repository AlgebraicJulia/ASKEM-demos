using Catlab, Catlab.Theories, Catlab.CategoricalAlgebra, Catlab.Graphs, Catlab.CategoricalAlgebra.FinSets
using Catlab.Graphics
using Catlab.CategoricalAlgebra.Diagrams
using Catlab.CategoricalAlgebra.FreeDiagrams
using Catlab.ACSetInterface
import ..FinCats: FreeCatGraph, FinDomFunctor, collect_ob, collect_hom
using Catlab.Programs

using PartialFunctions
using MLStyle

import Base.parse

#############################
# Some utils                 BEGIN
#############################

concatMap(f,xs) = foldl(vcat, map(f, xs))

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
monos(X::ACSet, Y::ACSet) = homomorphism(X,Y, monic=true)

"""Ask: "does there exists a mono X ↪ Y ?" 
"""
exists_mono(X::ACSet,Y::ACSet)::Bool = is_homomorphic(X,Y, monic=true)


#############################
# Maximum Common Acset (MCA)    BEGIN
#############################
"""Brute-force implementation of Maximum Common Acset (MCA).
Input: two Acsets a₁ and a₂ and an integer k (which is assumed to be zero if unspecified)
Task : find all a with with |a| ≥ k such that there is a monic span of Acset a₁ ← a → a₂. 
"""
function mca(XX::ACSet, YY::ACSet, k=0)
  (X,Y) = smallFirst(XX,YY)
  exists_mono(X,Y) ? (size(X) ≥ k ? [X] : []) : mca_help(X, Y, k)
end 

function mca_help(X::ACSet, Y::ACSet, k)
  C = acset_schema(X) #X: C → Set
  getsmaller(c) = map(dom, one_removed_subobs(X, c))
  # enumerate all sub-acsets X' ↪ X of the acset X: C → Set obtained by removing one point from Xc for some c ∈ C
  oneRemovedFromX = filter(
                      s -> size(s) ≥ k, #only keep the big ones
                      concatMap(getsmaller, filter(c -> !isempty(parts(X,c)), objects(C)))
                    )
  #keep only those X' ↪ X such that there exists mono X' ↪ Y
  Y_subs          = filter(χ -> exists_mono(χ, Y), oneRemovedFromX)
  #either terminate or recurse
  isempty(Y_subs) ? concatMap(χ -> mca_help(χ, Y), Y_subs) : Y_subs
end
#############################
# Maximum Common Acset (MCA)    END
#############################