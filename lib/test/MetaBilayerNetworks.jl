using ASKEM
using ASKEM.MetaBilayerNetworks
using Catlab
using Catlab.CategoricalAlgebra
using Catlab.Graphics
using AlgebraicPetri
using AlgebraicPetri.BilayerNetworks

bln = @acset MetaBilayerNetwork begin
  Qin = 2
  Qout = 2
  Box = 1
  Win = 2
  Wa = 2
  Wn = 2

  arg = [1,2]
  call = [1,1]

  influx = [1,1]
  efflux = [1,1]

  infusion = [2,2]
  effusion = [1,2]

  parameter = [:β]
  variable = [:S, :I]
  tanvar = [:S, :I]

  bmeta = ["the infection reaction"]
  qimeta = [(onto="susceptible people",unit="ppm"), (onto="infected people (ppm)", unit="ppm")]
  qometa = [(onto="susceptible people flux",unit="ppm/s"), (onto="infected people flux", unit="ppm/s")]
end

to_graphviz(bln)