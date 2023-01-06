module MetaBilayerNetworks

export SchMetaBilayerNetwork, MetaBilayerNetwork, MetaBilayerNetworkUntyped, AbstractMetaBilayerNetwork

using Catlab
using Catlab.Syntax
using Catlab.Theories
using Catlab.CategoricalAlgebra
using Catlab.Graphics
using AlgebraicPetri
using AlgebraicPetri.BilayerNetworks

@present SchMetaBilayerNetwork <: ThLabelledBilayerNetwork begin
  Metadata::AttrType

  qimeta::Attr(Qin, Metadata)
  qometa::Attr(Qout, Metadata)
  bmeta::Attr(Box, Metadata)
end

@abstract_acset_type AbstractMetaBilayerNetwork <: AbstractLabelledBilayerNetwork
@acset_type MetaBilayerNetworkUntyped(SchMetaBilayerNetwork) <: AbstractMetaBilayerNetwork
const MetaBilayerNetwork = MetaBilayerNetworkUntyped{Symbol, Any}
end