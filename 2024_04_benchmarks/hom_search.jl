"""
This file showcases the improvements to homomorphism search recently added to 
AlgebraicJulia. Our example homomorphism search problem will be the subgraph
isomorphism problem, i.e. looking for all matches of a (relatively small) 
pattern graph that live within a (relatively large) host graph.
"""

using Catlab, BenchmarkTools


c3 = cycle_graph(Graph, 3)

prog1 = compile_hom_search(c3, strat=:neighbor);
prog2 = compile_hom_search(c3, strat=:connected);

for (i, G) in [(i, erdos_renyi(Graph, 60, 0.05)) for i in 1:5]
  @btime homomorphisms($c3, $G)
  @btime homomorphisms($c3, $G; alg=VMSearch(), prog=$prog1)
  @btime homomorphisms($c3, $G; alg=VMSearch(), prog=$prog2)
end

