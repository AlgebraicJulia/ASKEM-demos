module TestMira 

using Test 
using ASKEM

data = joinpath(@__DIR__,"../../data/")
 
mdl_orig = load_mira("$data/BIOMD0000000971_petri_orig.json")
mdl_orig = load_mira_curated("$data/BIOMD0000000971_petri_curated.json")


end 