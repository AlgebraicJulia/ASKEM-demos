#!/usr/bin/env julia
#
# ./validatebilayer.jl bilayer.json petri.json will convert a bilayer to a petri for validation purposes.
using Catlab
using AlgebraicPetri
using AlgebraicPetri.BilayerNetworks

inputfile = ARGS[1]
outputfile = ARGS[2]

bln = read_json_acset(BilayerNetwork, inputfile)
show(bln)
lrn = LabelledReactionNet()
migrate!(lrn, bln)
write_json_acset(lrn, outputfile)