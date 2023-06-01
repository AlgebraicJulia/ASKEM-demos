module SubACSets
export mca, visualize_overlap

using Catlab, Catlab.Theories, Catlab.CategoricalAlgebra, Catlab.Graphs, Catlab.CategoricalAlgebra.FinSets
using Catlab.Graphics
using Catlab.CategoricalAlgebra.Diagrams
using Catlab.CategoricalAlgebra.FreeDiagrams
using Catlab.ACSetInterface
import Catlab.CategoricalAlgebra.FinCats: FreeCatGraph, FinDomFunctor, collect_ob, collect_hom
using Catlab.Programs

using AlgebraicPetri
using AlgebraicPetri: Subgraph, Attributes, Edge, Graphviz
using ..Ontologies

concatmap(f,xs) = mapreduce(f, vcat, xs; init =[])

# the smaller acset is the "pattern"
function normalize_order(X::ACSet, Y::ACSet)
  size(X) ≤ size(Y) ? (X,Y) : (Y,X)
end

# given an acset X: C→Set and an object c ∈ C, compute all possible X' ↪ X which
# are given by removing a single element from the set Xc
function one_removed_subobs(X::ACSet, c)
  C    = acset_schema(X)
  subs = [filter(j -> j != i, parts(X,c)) for i ∈ parts(X,c)]
  build_sub(sub) = hom(
    # need double negation to remove dangling edges
    ¬¬Subobject(
      X,
      NamedTuple(
        zip(
          objects(C),
          [o == c ? sub : parts(X, o) for o ∈ objects(C)]
        )
      )
    )
  )
  map(build_sub, subs)
end

"""
Defintion: let 𝐺: C → 𝐒et be a C-set, we define the _size_ of 𝐺 to be ∑_{c ∈ C}
|𝐺c|.  For example, under this definition, the size of:
  * a graph G is |GE| + |GV| (num edges + num vertices)
  * a Petri net P is |PT| + |PS| + |PI| + |PO| (num transitions + num species +
    num input arcs + num output arcs).
"""
function size(X::ACSet)
  foldl(+, [length(parts(X, oₛ)) for oₛ ∈ objects(acset_schema(X))])
end


"""Get all monomorphisms from an acset X to an acset Y
"""
monos(X::ACSet, Y::ACSet) =
  homomorphism(X, strip_names(Y); monic = true, type_components=(Name=x->nothing,),)

"""Ask: "does there exists a mono X ↪ Y ?"
"""
exists_mono(X::ACSet,Y::ACSet)::Bool =
  is_homomorphic(X, strip_names(Y); monic = true, type_components=(Name=x->nothing,),)


"""Brute-force implementation of Maximum Common Acset (MCA).
Input: two Acsets a₁ and a₂
Task : find all a with with |a| maximum possible such that there is a monic span of Acset a₁ ← a → a₂.
"""
function mca(XX::ACSet, YY::ACSet)
  (X,Y) = normalize_order(XX,YY)
  exists_mono(X,Y) ? [X] : mca_help(X, Y)
end

function mca_help(X::ACSet, Y::ACSet)
  C = acset_schema(X) #X: C → Set
  getsmaller(c) = map(dom, one_removed_subobs(X, c))
  # enumerate all sub-acsets X' ↪ X of the acset X: C → Set obtained by removing one point from Xc for some c ∈ C
  oneRemovedFromX = concatmap(getsmaller, filter(c -> !isempty(parts(X,c)), objects(C)))
  # keep only those X' ↪ X such that there exists mono X' ↪ Y
  Y_subs          = filter(χ -> exists_mono(χ, Y), oneRemovedFromX)
  # either terminate or recurse
  if isempty(Y_subs)
    concatmap(χ -> mca_help(χ, Y), Y_subs)
  else
    ω = maximum(map(size,Y_subs))
    filter(y -> size(y) == ω, Y_subs)
  end
end


"""
    O --> B
    ↓     ↓
    A -->⌜A+B         

Visualize the relation between two models via drawing the induced cospan of the 
pushout.
"""
visualize_overlap(o_span) = o_span |> colimit |> cocone |> Graph 

function Graph(m::Multicospan{<:AbstractPetriNet})
  apex = Subgraph(m.apex; name="clusterApex", pre="ap_")

  leg_graphs = map(enumerate(m.legs)) do (i, l)
    Subgraph(l.dom; name="clusterLeg$i", pre="leg$(i)_")
  end

  graph_attrs = Attributes(:rank=>"same")
  node_attrs  = Attributes(:shape=>"plain", :style=>"filled", :color=>"white")
  edge_attrs  = Attributes(:splines=>"splines")

  t_edges = map(enumerate(m.legs)) do (i, l)
    map(1:nt(dom(l))) do t
      Edge(["\"leg$(i)_t$(t)\"", "\"ap_t$(l.components[:T](t))\"",],
      Attributes(:constraint=>"false", :style=>"dotted"))
    end
  end

  s_edges = map(enumerate(m.legs)) do (i, l)
    map(1:ns(dom(l))) do s
      Edge(["\"leg$(i)_s$(s)\"", "\"ap_s$(l.components[:S](s))\"", ],
      Attributes(:constraint=>"false", :style=>"dotted"))
    end
  end
  if length(leg_graphs) == 2
    stmts = vcat([first(leg_graphs),apex,last(leg_graphs)], t_edges..., s_edges...)
  else 
    stmts = vcat([apex], leg_graphs, t_edges..., s_edges...)
  end 
  g = Graphviz.Digraph("G", stmts; graph_attrs=graph_attrs, node_attrs=node_attrs, edge_attrs=edge_attrs)
end

end
