# Benchmarks April 2024

AlgebraicJulia is structured around a small number of core categorical concepts and data structures, thus performance improvements to this core has a large multiplying force throughout the ecosystem. 

## VM Homomorphism search

Morphism search is core to many algorithms, such as finding possible pattern matches for applying rewrite rules. 

Our example homomorphism search problem will be the subgraph isomorphism problem, i.e. looking for all matches of a (relatively small) 
pattern graph that live within a (relatively large) host graph.

The algorithmic improvement over the default method of homomorphism search ("backtracking search") is called the Virtual Machine approach. In this case, before searching for homomorphisms $X\rightarrow Y$, we precompile a program optimized purely for finding morphisms out of $X$, which has various performance advantages, such as being able to allocate a fixed amount of memory for the search problem before running the algorithm (in contrast to dynamically allocating during the search process). 

The following data comes from looking for 3-cycles in a 60-vertex random graph. (Times given in ms)

| Run | Backtracking | VM1 | VM2 | Speedup |
|---|---|---|---|---|
| 1 | 257 | 28 | 22 | 12x | 
| 2 | 266 | 28 | 22 | 12x | 
| 3 | 263 | 31 | 25 | 11x |
| 4 | 355 | 38 | 31 | 11x | 
| 5 | 218 | 24 | 19 | 11x |

## In-place Rewriting

Our general framework for rewriting uses the general categorical concept of a pushout. This operaiton is most naturally specified in a pure, functional style, which does not mutate the original data. However, we have implemented an in-place, mutating version.

The following benchmark performs a rewrite rule on weighted graphs which replaces a pair of parallel edges with a single edge weighted by the sum, applied to a random graph of 50 vertices. Times are given in ms.

| Run | Functional | Inplace | Speedup |
|---|---|---|---|
| 1 | 4814 | 101.5| 47x |
| 2 | 4915 | 99.8 | 49x | 
| 3 | 4907 | 97.7 | 50x |

## Incremental Hom Search

Simulations require randomly sampling from sets of possible actions, and a key computational challenge is maintaining the set of possible actions as the world is changing. Although one can look at the state of the world and compute the possible actions from scratch (an operation which scales with the size of the simulated world, which is usually big), it is more efficient to just analyze _how the world changed_ and then incrementally update the current set of possible actions. The latter option scales with just the size of the change made to the world, which is small. The methodology is described in [this blog post](https://www.localcharts.org/t/incremental-presheaf-hom-set-updating/13224). 

This benchmark compares the start-from-scratch approach to the incremental update approach in a 200 vertex graph, after applying a random rewrite rule (involving appx 3-5 vertex graphs): 

| Run | From scratch | Incremental | Speedup |
|---|---|---|---|
| 1 | 1764 | 1.7 | 1037x |
| 2 | 7472 | 1.3 | 5748x |
| 3 | 4014| 10.2 | 394x |

(Times above in ms)

In general any time the world state is large there should be massive speed gains.

## Overall agent-based model workflow

A 100-step ABM simulation with 140 sheep and 20 wolves on a 30x30 grid took on average **3.012 s**. In Agents.jl, this took on average **0.94 ms**. Agents.jl is a purely imperative framework for modeling, in contrast to AlgebraicRewriting's purely declarative syntax. 

Caveat: this only leverages the in-place rewriting, not the VM nor incremental hom search improvements. 
