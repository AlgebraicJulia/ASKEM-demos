
# Rewriting in-place
using AlgebraicRewriting, Catlab
# Sum parallel edges of Weighted Graph 

# Datatypes
@acset_type MADWeightedGraph(SchWeightedGraph, part_type=BitSetParts, 
                             index=[:src,:tgt]) <: AbstractGraph
@acset_type WeightedGraph′(SchWeightedGraph, part_type=DenseParts, 
                             index=[:src,:tgt]) <: AbstractGraph

const MADIntGraph = MADWeightedGraph{Int}
const IntGraph = WeightedGraph′{Int}

# Specific graphs
L = @acset MADIntGraph begin 
  V=2; E=2; Weight=2; src=1; tgt=2; weight=AttrVar.(1:2) 
end

R = @acset MADIntGraph begin 
  V=2; E=1; Weight=1; src=1; tgt=2; weight=[AttrVar(1)] 
end

L′, R′ = IntGraph(), IntGraph();
copy_parts!.([L′, R′], [L, R])

# Rule
l, r = homomorphism.(Ref(MADIntGraph(2)), [L, R]; monic=true)
l′, r′ = homomorphism.(Ref(IntGraph(2)), [L′, R′]; monic=true)

plus(xs) = xs[1] + xs[2]

rule, rule′ = Rule.([l,l′], [r,r′]; monic=[:E], expr=Dict(:Weight=>[plus]))

# Apply rewrite
prog = compile_rewrite(rule)

for _ in 1:3
  G = erdos_renyi(MADIntGraph, 50, .3)
  G[:weight] = rand(1:10, nparts(G, :E))
  G′ = IntGraph();
  copy_parts!(G′, G)

  m, m′ = homomorphism.([L, L′], [G, G′])

  @btime rewrite_match!(rule, m; prog);
  @btime rewrite_match(rule′, m′);
end

