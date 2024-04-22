"""Incremental hom search benchmarks"""


using Catlab, AlgebraicRewriting
using AlgebraicRewriting.Rewrite.Utils: get_result, get_rmap, get_pmap
using AlgebraicRewriting.Incremental: validate, connected_acset_components, 
                                      state, deletion!, addition!
all_graph(rng) = [cycle_graph.(Graph,rng); star_graph.(Graph, rng);
                  path_graph.(Graph, rng)]
rand_graph(rng) = rand(all_graph(rng))

while true
  L, R = rand_graph(3:5), rand_graph(3:5)
  I = rand(path_graph.(Graph, 2:4))
  NV=200
  start = erdos_renyi(Graph, NV, 2*NV)
  l = homomorphism(I, L; monic=true)
  r = homomorphism(I, R; monic=true)
  isnothing(r) && continue
  m = homomorphism(I, start)
  isnothing(m) && continue
  res = rewrite_match_maps(Rule(id(I), r), m);
  (pl, pr), rmap = get_pmap(:DPO, res), get_rmap(:DPO, res);
  @assert collect(pl[:V]) == 1:NV

  @btime begin 
    new_matches = homomorphisms(L, codom($rmap))
  end;
  hset = IncHomSet(L, [r], start);
  @btime begin 
    deletion!($hset, $pl);
    addition!($hset, 1, $rmap, $pr);
  end 
  break
end