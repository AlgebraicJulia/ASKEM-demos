# module TestSubACSets

# using Test

using AlgebraicPetri
using Catlab.CategoricalAlgebra, Catlab.Graphics

# Kris Brown's Catlab commit: fd4427c04b4cbb00ef8d4d2a8bf19a05fb8b3eac

X = LabelledPetriNet(
  [:X, :Ya, :Wa, :Z],
  :f => (:X => (:Wa, :Ya, :Z))
)

Y = LabelledPetriNet(
  [:Wb, :X, :Yb, :Z],
  :f => ((:Wb, :X, :Yb) => :Z)
)

X = LabelledPetriNet(
  [:X, :Y, :W, :Z],
  :f => (:X => (:W, :Y, :Z))
)

Y = LabelledPetriNet(
  [:W, :X, :Y, :Z],
  :f => ((:W, :X, :Y) => :Z)
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
    # push!(res,¬¬Subobject(X,NamedTuple(zip([:T,:S],[X_sub_t,X_sub_s]))))
    push!(res,Subobject(X,NamedTuple(zip([:T,:S],[X_sub_t,X_sub_s]))))
  end 

  res
end
function mca_jank_af(X_list)
  As = intersect([X[:sname] for X in X_list]...)
  At = intersect([X[:tname] for X in X_list]...)
  # C = acset_schema(X) #X: C → Set

  res = []
  for X in X_list
    X_sub_s = [incident(X,s,:sname)[1] for s in As] 
    X_sub_t = [incident(X,t,:tname)[1] for t in At]

    # push!(res,¬¬Subobject(X,NamedTuple(zip([:T,:S],[X_sub_t,X_sub_s]))))
    push!(res,Subobject(X,NamedTuple(zip([:T,:S],[X_sub_t,X_sub_s]))))
  end 

  res
end

function mca_jank_af(X_list)
  As = intersect([X[:sname] for X in X_list]...)
  At = intersect([X[:tname] for X in X_list]...)
  # C = acset_schema(X) #X: C → Set

  state_props = Dict(Symbol(s["id"]) => s for s in model["states"])
  transition_props = Dict(Symbol(t["id"]) => t for t in model["transitions"])
  transitions = [Symbol(t["id"]) => (Symbol.(t["input"]) => Symbol.(t["output"])) for t in model["transitions"]]

  # A = PropertyLabelledPetriNet{Dict}(LabelledPetriNet(As, transitions...), state_props, transition_props)
  A = LabelledPetriNet()
  add_parts!(A,:S,length(As);sname=As)
  add_parts!(A,:T,length(At);tname=At)
  Ax = []
  for X in X_list
    S = [incident(X,s,:sname)[1] for s in As]
    T = [incident(X,t,:tname)[1] for t in At]
    I = []
    O = []
    push!(Ax,LooseACSetTransformation((S=S, T=T, I=I, O=O), (Name=x->nothing), A, X))

  end

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

test = pushout(hom(A[1]),hom(A[2]))
intersect(legs(test))[1] |> to_graphviz

#=
Should directly take intersection of PetriNets
Want to take intersection of double negations
so need to be of subobjects of common acset_schema
so take pushout of inclusion of discrete A into the double negations

mca(P,Q,sid,tid)
sid: AbstractPetri --> Vector{T} (e.g. T is string or symbol)
tid: AbstractPetri --> Vector{T}

default: sid(P) --> P[:,:sname]
o.w. id property 
=#
# end

m1 = LabelledPetriNet(
  [:X1, :Y1, :W1, :Z1],
  :f1 => (:X1 => (:W1, :Y1, :Z1))
)

m2 = LabelledPetriNet(
  [:X2, :W2, :Y2, :Z2],
  :f2 => ((:W2, :X2, :Y2) => :Z2)
)

mca1, mca1_morphs = mca(m1, m2)

@test is_isomorphic(mca1,strip_attributes(LabelledPetriNet(
  [:X1, :Y1, :W1, :Z1],
  :f1 => (:X1 => :Z1)
)))
@test length(mca1_morphs) == 2
@test length(mca1_morphs[1]) == 3
@test length(mca1_morphs[2]) == 3


@test is_isomorphic(PetriNet(mca1),mca(PetriNet(m1), PetriNet(m2))[1])

m3 = LabelledPetriNet(
  [:W3, :X3, :Y3, :Z3],
  :f3 => ((:W3, :X3) => (:Y3, :Z3))
)

mca3, mca3_morphs = mca([m3, m2, m1])
@test is_isomorphic(mca3,strip_attributes(LabelledPetriNet(
  [:W3, :X3, :Y3, :Z3],
  :f3 => (:X3 => :Z3)
)))
@test length(mca3_morphs) == 3
@test length(mca3_morphs[1]) == 4
@test length(mca3_morphs[2]) == 3
@test length(mca3_morphs[3]) == 3

m4 = LabelledPetriNet(
  [:X4, :Y4],
  :f41 => (:X4 => :Y4),
  :f42 => (:Y4 => :Y4)  
)
mca4, mca4_morphs = mca([m3, m2, m1, m4])
@test is_isomorphic(mca4,strip_attributes(LabelledPetriNet(
  [:X4, :Y4],
  :f41 => (:X4 => :Y4)
)))
@test length(mca4_morphs) == 4
#=
@test length(mca4_morphs[1]) == 4
@test length(mca4_morphs[2]) == 3
@test length(mca4_morphs[3]) == 3
@test length(mca4_morphs[4]) == 2
=#

m5 = LabelledPetriNet(
  [:X5, :Y5],
  :f51 => (:X5 => :Y5),
  :f52 => (:Y5 => :X5)  
)
mca5, mca5_morphs = mca([m3, m2, m1, m5])

using OrderedCollections
import Catlab.CategoricalAlgebra.CSets: maximum_common_subobject, total_parts
function maximum_common_subobject(Xs::Vector{T}) where T <: ACSet
  it = partial_overlaps(Xs)
  print(length(collect(it)))
  osize = -1
  res = OrderedDict()
  for overlap in it 
    apx = apex(overlap)
    size = total_parts(apx)
    osize = osize == -1 ? size : osize
    println(size,", ",osize)
    if size < osize return res end 
    res[apx] = overlap
  end 
  return res
end