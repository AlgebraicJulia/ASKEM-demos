using Catlab
using Catlab.Syntax
using Catlab.Theories
import JSON, JSONSchema

using Catlab, Catlab.Theories, Catlab.Graphs, Catlab.CategoricalAlgebra
using Catlab.Graphics

@present SchMira(FreeSchema) begin
    (V, P, D, NC, CC, GCC)::Ob
    # attributes for rate number on each template element
    # rate law XML blob 

    subject_P::Hom(P, V)
    outcome_D::Hom(D, V)

    subject_NC::Hom(NC, V)
    outcome_NC::Hom(NC, V)

    subject_CC::Hom(CC, V)
    control_CC::Hom(CC, V)
    outcome_CC::Hom(CC, V)

    subject_GCC::Hom(GCC, V)
    control₁_GCC::Hom(GCC, V)
    control₂_GCC::Hom(GCC, V)
    outcome_GCC::Hom(GCC, V)
end

@present SchNamedMira <: SchMira begin
    Name::AttrType
    nameV::Attr(V, Name)
    nameP::Attr(P, Name)
    nameD::Attr(D, Name)
    nameNC::Attr(NC, Name)
    nameCC::Attr(CC, Name)
    nameGCC::Attr(GCC, Name)
end

@present SchMiraModel <: SchNamedMira begin
    (Rate, Ex)::AttrType
    rateP::Attr(P, Rate)
    rateD::Attr(D, Rate)
    rateNC::Attr(NC, Rate)
    rateCC::Attr(CC, Rate)
    rateGCC::Attr(GCC, Rate)

    lawP::Attr(P, Ex)
    lawD::Attr(D, Ex)
    lawNC::Attr(NC, Ex)
    lawCC::Attr(CC, Ex)
    lawGCC::Attr(GCC, Ex)
end

@acset_type MiraModel(SchMiraModel, index=[])

M = @acset MiraModel{Symbol, Float64, Symbol} begin
    V = 4
    CC = 1 
    NC = 1
    P = 1
    D = 1
    GCC=1

    subject_P = [1]
    outcome_D = [1]
    subject_NC = [2]
    outcome_NC = [3]
    subject_CC = [2]
    control_CC = [3]
    outcome_CC = [1]

    subject_GCC = [2]
    control₁_GCC = [3]
    control₂_GCC = [4]
    outcome_GCC = [1]

    nameV = [:v1, :v2, :v3, :v4]
    nameP = :p
    nameD = :d
    nameNC = :nc
    nameCC = :cc
    nameGCC = :gcc

    rateP =  1.0
    rateD =  1.0
    rateNC = 1.0
    rateCC = 1.0
    rateGCC = 1.0

    lawD = :const
    lawP = :const
    lawNC = :linear
    lawCC = :MAK
    lawGCC = :generic
end

add_part!(M, :NC, subject_NC=4, outcome_NC=1, nameNC=:return, rateNC=2.0, lawNC=:linear)
M
JSON.print(M, 2)