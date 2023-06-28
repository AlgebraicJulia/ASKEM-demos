# module TestSubACSets

# using Test

using AlgebraicPetri
using Catlab.CategoricalAlgebra, Catlab.Graphics

X = LabelledPetriNet(
  [:X, :Ya, :Wa, :Z],
  :f => (:X => (:Wa, :Ya, :Z))
)

Y = LabelledPetriNet(
  [:Wb, :X, :Yb, :Z],
  :f => ((:Wb, :X, :Yb) => :Z)
)

function mca_jank_af(X,Y)
  mca_jank_af([X,Y])
end

function mca_jank_af(X_list)
  As = intersect([X[:sname] for X in X_list]...)
  At = intersect([X[:tname] for X in X_list]...)
  # C = acset_schema(X) #X: C → Set

  res = []
  for X in X_list
    X_sub_s = filter(j -> X[:sname][j] ∈ As, parts(X,:S)) 
    X_sub_t = filter(j -> X[:tname][j] ∈ At, parts(X,:T))

    # hom(¬¬Subobject(X,NamedTuple(zip(objects(C),[o == c ? sub : parts(X, o) for o ∈ objects(C)]))))
    # hom(¬¬Subobject(X,NamedTuple(zip(objects(C),[X_sub_t,X_sub_s,parts(X,:I),parts(X,:O)]))))
    push!(res,¬¬Subobject(X,NamedTuple(zip([:T,:S],[X_sub_t,X_sub_s]))))
  end

  res
end

A = mca_jank_af(X,Y)
A[1] |> to_graphviz
A[2] |> to_graphviz

dom(hom(A[1]))==dom(hom(A[2]))

# end
