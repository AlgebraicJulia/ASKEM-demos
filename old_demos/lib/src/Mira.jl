module Mira 
export load_mira, load_mira_curated

using JSON
using AlgebraicPetri
using Catlab, Catlab.CategoricalAlgebra

@present SchOrigMIRANet <: SchLabelledReactionNet begin
  MID::AttrType
  MCTX::AttrType
  Template::AttrType
  mira_ids::Attr(S, MID)
  mira_context::Attr(S, MCTX)
  template_type::Attr(T, Template)
  parameter_name::Attr(T, Name)
  parameter_value::Attr(T, Rate)
end
@abstract_acset_type AbstractOrigMIRANet <: AbstractLabelledReactionNet
@acset_type OrigMIRANet(SchOrigMIRANet) <: AbstractOrigMIRANet

@present SchMIRANet <: SchLabelledReactionNet begin
  MID::AttrType
  MCTX::AttrType
  Template::AttrType
  mira_ids::Attr(S, MID)
  mira_context::Attr(S, MCTX)
  template_type::Attr(T, Template)
  parameter_name::Attr(T, Name)
  # parameter_value::Attr(T, Rate)
end
@abstract_acset_type AbstractMIRANet <: AbstractLabelledReactionNet
@acset_type MIRANet(SchMIRANet) <: AbstractMIRANet


function load_mira(fname)
  mdl_orig_mira = read_json_acset(OrigMIRANet{Any,Any,Any,Any,Any,Any}, fname)
  mdl_orig_mira[:sname] .= Symbol.(mdl_orig_mira[:sname])
  mdl_orig_mira[:tname] .= Symbol.(mdl_orig_mira[:tname])
  mdl_orig = LabelledPetriNet(mdl_orig_mira);
  return mdl_orig
end

function load_mira_curated(fname, fillrate=1.0)
  mdl_disease_mira = read_json_acset(MIRANet{Any,Any,Any,Any,Any,Any}, fname)
  input_rates = deepcopy(mdl_disease_mira[:rate])
	min_rate = minimum(input_rates[map(!isnothing,input_rates)]) # UNUSED!
  filled_rates = input_rates
  filled_rates[map(isnothing,input_rates)] = repeat([fillrate],sum(map(isnothing,input_rates)))
  mdl_disease_mira[:rate] = filled_rates
  mdl_disease_mira[:concentration] = 0.0
  mdl_disease_mira[:sname] .= Symbol.(mdl_disease_mira[:sname])
  mdl_disease_mira[:tname] .= Symbol.(mdl_disease_mira[:tname])
  return mdl_disease_mira
end


end # module 
