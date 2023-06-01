#!/usr/bin/env julia
#
# ./validatebilayer.jl bilayer.json petri.json will convert a bilayer to a petri for validation purposes.
#
# Note that, to run, the --project switch must point to a directory 
# containing a Project.toml file with Catlab and AlgebraicPetri added.
# For example, to run from the ASKEM-demos/scripts directory:
# julia --project="../Dec2022Demo" -- validatebilayer.jl ../data/CHIME_SVIIvR_dynamics_BiLayer.json chime_sviivr_lpn.json
#
using Catlab
using Catlab.CategoricalAlgebra
using AlgebraicPetri
using AlgebraicPetri.BilayerNetworks

inputfile = ARGS[1]
outputfile = ARGS[2]

bln = read_json_acset(LabelledBilayerNetwork, inputfile)
show(bln)
# lrn = LabelledReactionNet()
# migrate!(lrn, bln)
# write_json_acset(lrn, outputfile) 
lpn = LabelledPetriNet()
migrate!(lpn, bln)
write_json_acset(lpn, outputfile)